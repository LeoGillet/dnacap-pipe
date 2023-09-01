from src import amr
from src import mlst


def extra_options():
    return "\n".join([
        "\t11. Start AMR detection",
        "\t12. Start MLST typing"
        ])


def start(step: int):
    match step:
        case 11:
            amr.start_pylori_amr_analysis()
        
        case 12:
            mlst.start_mlst_typing('pylori')
