from src import amr
from src import mlst


def extra_options():
    return "".join([
        "\t11. Start AMR detection",
        "\t12. Start MLST typing"
        ])


def start(step: int):
    match step:
        case 11:
            muts = amr.start_pylori_amr_analysis()
        
        case 12:
            print(mlst._concatenate_pylori_mlst("DC2_10_S19"))
