"""Identidad química de los fármacos y anotación de conjugabilidad a Fe3O4@SiO2-APTES.

La identidad se resuelve contra PubChem por nombre; la correspondencia nombre->CID es
ambigua por construcción y se marca como tal. La anotación de conjugabilidad describe
grupos funcionales presentes, no predice actividad ni eficiencia de carga: la química
real se optimiza en el laboratorio para cada candidato.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = ("Title,MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,"
              "InChIKey,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,Charge")
# PubChem documenta un máximo de 5 peticiones por segundo; se usa un margen amplio.
DELAY_S = 0.25
TIMEOUT_S = 30
RETRIES = 3

# Patrones de grupos funcionales. Cada clave se anota como presencia/ausencia observada
# en la estructura; la accesibilidad estérica real no se evalúa aquí.
SMARTS = {
    # Ruta covalente principal: carboxilo libre + amina del APTES con EDC/NHS.
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    # Rutas covalentes alternativas: requieren activar la superficie (p. ej. glutaraldehído).
    "primary_aliphatic_amine": "[NX3;H2;!$(NC=[O,N,S]);!$(N[a])]",
    "secondary_aliphatic_amine": "[NX3;H1;!$(NC=[O,N,S]);!$(N[a])]",
    # Las aminas unidas a anillo aromático, y sobre todo a heterociclos tipo aminopirimidina,
    # son poco nucleófilas: no se acoplan en condiciones EDC/NHS habituales.
    "aromatic_amine": "[NX3;H2,H1][a]",
    # Asas derivatizables: requieren un paso previo (p. ej. succinilación) que genera el COOH.
    "aliphatic_hydroxyl": "[OX2H][CX4]",
    "phenol": "[OX2H][c]",
    "thiol": "[SX2H]",
    "boronic_acid": "[BX3]([OX2H])[OX2H]",
    # Grupos aniónicos a pH fisiológico: favorecen adsorción sobre superficie protonada.
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H1]",
    "phosphate_acid": "[PX4](=O)([OX2H1])",
    "tetrazole": "c1nnn[nH]1",
    "acyl_sulfonamide": "[CX3](=O)[NX3H][SX4](=O)(=O)",
    # Grupos básicos: a pH fisiológico quedan catiónicos, como la propia superficie.
    "amidine": "[NX3][CX3]=[NX2]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "basic_tertiary_amine": "[NX3;H0;!$(NC=[O,N,S]);!$(N[a]);!$(N=*)]",
    # Pasivos de estabilidad acuosa: hidrólisis durante carga, almacenamiento o liberación.
    "ester": "[CX3](=O)[OX2H0][#6;!$([CX3]=[OX1])]",
    "lactone": "[CX3;R](=O)[OX2;R][#6;R]",
    "anhydride": "[CX3](=O)[OX2H0][CX3](=O)",
    "beta_lactam": "C1(=O)NCC1",
    "epoxide": "[OX2;r3]1[#6;r3][#6;r3]1",
    "aldehyde": "[CX3H1](=O)[#6]",
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    "hydrazone_imine": "[CX3]=[NX2][NX3,OX2,#6]",
}
# Un carboxilo dentro de estos entornos no es un asa de conjugación libre.
ACID_EXCLUSIONS = ("[CX3](=O)[OX2H1][#6]",)
HYDROLYSIS_FLAGS = ("ester", "lactone", "anhydride", "beta_lactam", "epoxide",
                    "acyl_halide", "hydrazone_imine")
ANIONIC_FLAGS = ("carboxylic_acid", "sulfonic_acid", "phosphate_acid", "tetrazole",
                 "acyl_sulfonamide")
CATIONIC_FLAGS = ("primary_aliphatic_amine", "secondary_aliphatic_amine", "amidine",
                  "guanidine", "basic_tertiary_amine")
AMINE_FLAGS = ("primary_aliphatic_amine", "secondary_aliphatic_amine")
DERIVATIZABLE_FLAGS = ("aliphatic_hydroxyl", "phenol", "thiol", "boronic_acid")
# Metales con química de coordinación propia: el precedente de oxaliplatino pertenece aquí.
COORDINATION_METALS = {"Pt", "Ru", "Pd", "Au", "Gd", "As"}
HIGH_LOGP = 5.0


def _request(url):
    """Devuelve el JSON de PubChem, None si el nombre no existe, y falla si la red falla."""
    last = None
    for attempt in range(RETRIES):
        try:
            response = requests.get(url, timeout=TIMEOUT_S)
        except requests.RequestException as error:  # red intermitente, no ausencia del dato
            last = error
            time.sleep(1.5 * (attempt + 1))
            continue
        if response.status_code == 404:
            return None
        if response.status_code == 200:
            return response.json()
        last = RuntimeError(f"HTTP {response.status_code} en {url}")
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"PubChem no respondió tras {RETRIES} intentos: {last}")


def name_variants(name):
    """Formas alternativas del nombre, en orden de confianza decreciente.

    Las anotaciones de concentración o estereoquímica entre paréntesis impiden la
    búsqueda; se prueba el nombre limpio, pero se registra cuál funcionó para revisar.
    """
    variants = [name]
    without_parentheses = re.sub(r"\s*\([^)]*\)", "", name).strip()
    if without_parentheses and without_parentheses != name:
        variants.append(without_parentheses)
    hyphen_as_space = without_parentheses.replace("-", " ").strip()
    if hyphen_as_space and hyphen_as_space not in variants:
        variants.append(hyphen_as_space)
    return variants


def fetch_compound(name, cache_dir):
    """Resuelve un nombre contra PubChem con caché en disco; conserva la respuesta cruda."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)[:120]
    cached = cache_dir / f"{safe}.json"
    if cached.exists():
        return json.loads(cached.read_text())

    payload = {"query_name": name, "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    payload["status"] = "not_found"
    for variant in name_variants(name):
        cids = _request(f"{PUBCHEM}/name/{quote(variant)}/cids/JSON")
        time.sleep(DELAY_S)
        if cids is None:
            continue
        found = cids["IdentifierList"]["CID"]
        payload.update(n_cids=len(found), query_used=variant)
        properties = _request(f"{PUBCHEM}/cid/{found[0]}/property/{PROPERTIES}/JSON")
        time.sleep(DELAY_S)
        payload["status"] = "resolved" if properties else "cid_without_properties"
        if properties:
            payload["properties"] = properties["PropertyTable"]["Properties"][0]
        break
    cached.write_text(json.dumps(payload, indent=2))
    return payload


def largest_fragment(smiles):
    """Componente mayor de una sal o mezcla; es donde se evalúan los grupos funcionales."""
    RDLogger.DisableLog("rdApp.*")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, None
    fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=False)
    if len(fragments) == 1:
        return molecule, molecule
    largest = max(fragments, key=lambda f: f.GetNumHeavyAtoms())
    # Los fragmentos salen sin sanear para tolerar valencias atípicas de complejos
    # metálicos; sin sanear falta la información de anillos que usan los SMARTS.
    try:
        Chem.SanitizeMol(largest)
    except (Chem.AtomValenceException, Chem.KekulizeException, ValueError):
        Chem.FastFindRings(largest)
    return molecule, largest


def functional_groups(molecule):
    """Conteo de coincidencias por patrón sobre el fragmento principal."""
    counts = {}
    for label, pattern in SMARTS.items():
        query = Chem.MolFromSmarts(pattern)
        counts[label] = len(molecule.GetSubstructMatches(query)) if query is not None else 0
    free_acid = Chem.MolFromSmarts(SMARTS["carboxylic_acid"])
    matches = set(molecule.GetSubstructMatches(free_acid))
    for pattern in ACID_EXCLUSIONS:
        excluded = Chem.MolFromSmarts(pattern)
        if excluded is not None:
            for match in molecule.GetSubstructMatches(excluded):
                matches.discard(tuple(match[:3]))
    counts["free_carboxylic_acid"] = len(matches)
    return counts


def annotate(smiles):
    """Anota una estructura: grupos, carga aproximada a pH 7,4, pasivos y ruta propuesta.

    La carga es una aproximación por grupos ionizables, no un cálculo de pKa.
    """
    whole, fragment = largest_fragment(smiles)
    if whole is None or fragment is None:
        return {"parse_status": "unparsed"}

    record = {"parse_status": "ok",
              "n_fragments": len(Chem.GetMolFrags(whole)),
              "heavy_atoms": fragment.GetNumHeavyAtoms()}
    record["chemical_form"] = "salt_or_mixture" if record["n_fragments"] > 1 else "single_component"
    elements = {atom.GetSymbol() for atom in whole.GetAtoms()}
    record["metal_complex"] = bool(elements & COORDINATION_METALS)
    record["metals"] = ",".join(sorted(elements & COORDINATION_METALS))
    try:
        record["mw_fragment"] = round(Descriptors.MolWt(fragment), 2)
    except Exception:  # estructuras con valencias atípicas, frecuentes en complejos metálicos
        record["mw_fragment"] = None

    record.update(functional_groups(fragment))
    record["anionic_groups"] = sum(record[k] for k in ANIONIC_FLAGS)
    record["cationic_groups"] = sum(record[k] for k in CATIONIC_FLAGS)
    record["net_charge_ph74"] = record["cationic_groups"] - record["anionic_groups"]
    record["hydrolysis_liability"] = sum(record[k] for k in HYDROLYSIS_FLAGS) > 0
    record["hydrolysis_groups"] = ",".join(k for k in HYDROLYSIS_FLAGS if record[k])
    record["alternative_amine_handle"] = sum(record[k] for k in AMINE_FLAGS) > 0
    record["derivatizable_handle"] = sum(record[k] for k in DERIVATIZABLE_FLAGS) > 0
    # La superficie APTES está protonada a pH fisiológico: un fármaco catiónico sin
    # grupos aniónicos se repele en lugar de adsorberse.
    record["cationic_repulsion_risk"] = (record["net_charge_ph74"] > 0
                                         and record["anionic_groups"] == 0)

    record["conjugation_route"] = route(record)
    # Conjugar el único grupo ionizable altera la forma en que el fármaco interactúa con
    # su diana. No es un veredicto de pérdida de actividad: requiere revisión química.
    record["single_handle_review"] = (record["free_carboxylic_acid"] == 1
                                      and record["anionic_groups"] == 1
                                      and record["cationic_groups"] == 0)
    return record


def route(record):
    """Ruta de unión propuesta según los grupos presentes. Prioriza covalente por EDC/NHS.

    El orden refleja cuánto trabajo químico adicional exige cada opción, no eficacia.
    """
    if record["free_carboxylic_acid"] > 0:
        return "covalente_cooh_edc_nhs"
    if record["metal_complex"]:
        return "coordinacion_metal"
    if record["anionic_groups"] > 0:
        return "adsorcion_electrostatica"
    if record["alternative_amine_handle"]:
        return "covalente_amina_requiere_activacion"
    if record["derivatizable_handle"]:
        return "covalente_requiere_derivatizacion"
    if record["aromatic_amine"] > 0:
        return "asa_aromatica_poco_reactiva"
    return "sin_asa_reactiva"


def build_table(names, cache_dir):
    """Tabla por nombre consultado con identidad, propiedades y anotación estructural."""
    rows = []
    for name in names:
        payload = fetch_compound(name, cache_dir)
        row = {"drug_name": name, "resolution_status": payload["status"],
               "n_pubchem_cids": payload.get("n_cids"), "source": "PubChem PUG REST",
               "retrieved_utc": payload["retrieved_utc"]}
        properties = payload.get("properties", {})
        row["pubchem_cid"] = properties.get("CID")
        row["pubchem_title"] = properties.get("Title")
        row["canonical_smiles"] = properties.get("SMILES") or properties.get("ConnectivitySMILES")
        row["inchikey"] = properties.get("InChIKey")
        row["molecular_formula"] = properties.get("MolecularFormula")
        row["molecular_weight"] = properties.get("MolecularWeight")
        row["xlogp"] = properties.get("XLogP")
        row["tpsa"] = properties.get("TPSA")
        row["hbd"] = properties.get("HBondDonorCount")
        row["hba"] = properties.get("HBondAcceptorCount")
        row["formal_charge"] = properties.get("Charge")
        # El nombre consultado y el título devuelto deben revisarse manualmente si difieren.
        row["name_matches_title"] = (isinstance(row["pubchem_title"], str)
                                     and row["pubchem_title"].strip().lower() == name.strip().lower())
        row["ambiguous_name"] = bool(row["n_pubchem_cids"] and row["n_pubchem_cids"] > 1)
        row["query_used"] = payload.get("query_used")
        row["query_was_modified"] = bool(row["query_used"] and row["query_used"] != name)
        if row["canonical_smiles"]:
            row.update(annotate(row["canonical_smiles"]))
            row["high_logp"] = (row["xlogp"] is not None and row["xlogp"] >= HIGH_LOGP)
        else:
            row["parse_status"] = "no_structure"
        rows.append(row)
    table = pd.DataFrame(rows)
    # Identidad química compartida: el bloque de conectividad del InChIKey agrupa
    # sales y alias del mismo esqueleto molecular.
    table["compound_id"] = table["inchikey"].str.slice(0, 14)
    return table
