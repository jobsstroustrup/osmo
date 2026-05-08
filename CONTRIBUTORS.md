# Contributors

## Core Team

- **Zhenyu He** ([@jobsstroustrup](https://github.com/jobsstroustrup))
  - Layer 1 (`hzy3.py`): designed and implemented the search mechanism, including relative coordinates, angular enumeration, torus periodic-boundary expansion, and 200-frame forward-simulation rollout verification
  - Layer 2 (`mechanic.py`): designed and implemented the parameter-tuning variant (sole author), exploring time-unit alignment, acceleration discretization, and temporal decay factor sensitivity
  - Layer 3 (`strategy.py`): designed and implemented the high-level strategy orchestration (sole author), integrating reactive hazard avoidance with deliberative goal planning, density-adaptive search radius, and stateful execution with mid-plan abort
  - Layer 4 reverse-engineering: complete reverse-engineering of `kernel.py` physics (eject, absorb, collision resolution, momentum conservation) into the `hzy3.py` rollout simulator helpers (`myeject`, `myabsorb`, `myupdate`)
  - Period: Spring 2019

- **szx**, **bs**, and other teammates
  - Layer 5 sample/ AI library co-development: opponent AI variants used as benchmark ladder (`player_ai szx.py`, `player_szx 1.0.py`, `player_szx 1.1.py`, etc.)
  - Collaborative discussion, learning, and verification across all layers

## Tournament Result

The team's submission placed **2nd in the SESSDSA 2019 Spring course tournament** (Peking University, Data Structures and Algorithms, Prof. Chen Bin).

## Acknowledgments

This project was completed as a team course assignment where all members participated in learning, discussion, and verification before dividing implementation labor by component.

For detailed methodology, architecture, and capability analysis, see the project README.
