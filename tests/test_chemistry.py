"""Pruebas de la anotación química y del cribado de candidatos. Sin red: sólo caché."""
import json
import sys
from pathlib import Path

import pytest

from gbm_drug_response.chemistry import annotate, build_table, fetch_compound, route

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

# Estructuras de referencia tomadas de PubChem; el valor esperado es la ruta química.
KNOWN = {
    # Dos carboxilos libres: caso directo de acoplamiento EDC/NHS con la amina del APTES.
    "Methotrexate": ("CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)NC(CCC(=O)O)C(=O)O",
                     "covalente_cooh_edc_nhs"),
    # Precedente del laboratorio: complejo de platino, química de coordinación propia.
    "Oxaliplatin": ("C1CCC(C(C1)[NH-])[NH-].C(=O)(C(=O)O)O.[Pt+2]", "coordinacion_metal"),
    # Amina alifática primaria del azúcar: ruta covalente previa activación de superficie.
    "Doxorubicin": ("CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)"
                    "(C(=O)CO)O)N)O", "covalente_amina_requiere_activacion"),
    # Sin asa reactiva: la carboxamida no se acopla y no hay grupo ionizable libre.
    "Temozolomide": ("CN1C(=O)N2C=NC(=C2N=N1)C(=O)N", "sin_asa_reactiva"),
    # Aminopirimidina: amina aromática poco nucleófila, no es asa práctica.
    "Dabrafenib": ("CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C(=CC=C3)NS(=O)(=O)"
                   "C4=C(C=CC=C4F)F)F", "asa_aromatica_poco_reactiva"),
}


@pytest.mark.parametrize("name", sorted(KNOWN))
def test_known_routes(name):
    smiles, expected = KNOWN[name]
    assert annotate(smiles)["conjugation_route"] == expected


def test_aromatic_amine_is_not_counted_as_aliphatic_handle():
    record = annotate(KNOWN["Dabrafenib"][0])
    assert record["aromatic_amine"] > 0
    assert not record["alternative_amine_handle"]


def test_ester_is_not_a_free_carboxylic_acid():
    # Acetato de metilo: el carboxilo está esterificado y no es un asa de conjugación.
    record = annotate("COC(C)=O")
    assert record["free_carboxylic_acid"] == 0
    assert record["hydrolysis_liability"]


def test_charge_and_repulsion_against_protonated_surface():
    anionic = annotate("OC(=O)CCCC(=O)O")       # ácido adípico
    cationic = annotate("CCCCN")                 # butilamina
    assert anionic["net_charge_ph74"] < 0 and not anionic["cationic_repulsion_risk"]
    assert cationic["net_charge_ph74"] > 0 and cationic["cationic_repulsion_risk"]


def test_largest_fragment_drives_annotation_of_a_salt():
    # Clorhidrato: el contraión no debe cambiar los grupos funcionales detectados.
    free = annotate("CCCCN")
    salt = annotate("CCCCN.Cl")
    assert salt["chemical_form"] == "salt_or_mixture"
    assert salt["primary_aliphatic_amine"] == free["primary_aliphatic_amine"]


def test_unparsed_structure_is_reported_not_guessed():
    assert annotate("no-es-smiles")["parse_status"] == "unparsed"


def test_route_prefers_carboxyl_over_other_handles():
    both = {"free_carboxylic_acid": 1, "metal_complex": False, "anionic_groups": 1,
            "alternative_amine_handle": True, "derivatizable_handle": True, "aromatic_amine": 1}
    assert route(both) == "covalente_cooh_edc_nhs"


def test_fetch_uses_cache_without_network(tmp_path, monkeypatch):
    payload = {"query_name": "Fake", "status": "resolved", "retrieved_utc": "2026-01-01T00:00:00Z"}
    (tmp_path / "Fake.json").write_text(json.dumps(payload))
    monkeypatch.setattr("gbm_drug_response.chemistry._request",
                        lambda url: pytest.fail("No debe consultarse la red con caché disponible"))
    assert fetch_compound("Fake", tmp_path) == payload


def test_build_table_groups_aliases_by_inchikey(tmp_path, monkeypatch):
    # Dos nombres del mismo esqueleto molecular deben compartir compound_id.
    responses = {
        "DrugA": {"CID": 1, "Title": "DrugA", "SMILES": "CCCCN",
                  "InChIKey": "HQABUPZFAYXKJW-UHFFFAOYSA-N"},
        "DrugA hydrochloride": {"CID": 2, "Title": "DrugA hydrochloride", "SMILES": "CCCCN.Cl",
                                "InChIKey": "HQABUPZFAYXKJW-UHFFFAOYSA-M"},
    }
    monkeypatch.setattr("gbm_drug_response.chemistry.fetch_compound",
                        lambda name, cache: {"status": "resolved", "n_cids": 1,
                                             "retrieved_utc": "2026-01-01T00:00:00Z",
                                             "properties": responses[name]})
    table = build_table(list(responses), tmp_path)
    assert table["compound_id"].nunique() == 1
    assert set(table["chemical_form"]) == {"single_component", "salt_or_mixture"}


def test_required_loading_scales_with_potency_and_mass():
    from rank_candidates import NP_PPM_MAX, required_loading_pct
    # A 25 ppm, un fármaco de PM 500 y IC50 1 µM necesita 2% p/p de carga.
    assert required_loading_pct(1.0, 500.0, NP_PPM_MAX) == pytest.approx(2.0)
    # Oxaliplatino libre en T98G (~242 µM, PM 397) exige una carga imposible.
    assert required_loading_pct(242.0, 397.0, NP_PPM_MAX) > 100
