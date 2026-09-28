from g_functions.g_functions import gfunctions
from g_functions.phase_diagram import build_phase_diagram

__all__ = ["gfunctions", "build_phase_diagram"]


if __name__ == "__main__":
    print(build_phase_diagram().head())
