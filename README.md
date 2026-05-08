# Osmo: Game AI with Multi-Layer Decision Architecture

An AI strategy for the Osmo absorption game (inspired by Hemisphere Games' Osmos).
Players control a cell that grows by absorbing smaller objects while avoiding larger threats.

## Project Overview

- **Status**: Complete, 2019 to 2020
- **Type**: Course final project (Peking University, Data Structures and Algorithms, SESSDSA 2019 Spring)
- **Team**: Zhenyu He (Layer 1 + Layer 2 + Layer 3 design and implementation), szx, bs, and other teammates (Layer 5 sample/ AI library co-development)
- **Codebase**: ~6,400 lines total

## Tournament Result

The team's submission placed **2nd in the SESSDSA 2019 Spring course tournament**.

## Core Game Mechanics

- Complete information (full world state visible to both AIs)
- Turn-based alternating moves
- Physics: momentum conservation, mass conservation, multi-body collisions
- Win condition: last cell standing (or most mass absorbed)

## AI Architecture (5 Layers)

### Layer 1: Search Mechanism (`hzy3.py`, 577 lines)

Author: Zhenyu He.

- Relative coordinate system + angular search
- Physics-informed pruning (limits 2D action space to 1D angle interval)
- Torus periodic boundary expansion (9 regions, 4 best candidates)
- Forward simulation verification (200-frame rollout to check plan viability)
- Custom node-based data structure for decision candidates

This layer combines physics intuition (relative coords, torus topology) with forward simulation (AlphaGo-style rollout) for decision verification.

### Layer 2: Parameter Tuning Variant (`mechanic.py`, 570 lines)

Author: Zhenyu He.

- Same search skeleton as Layer 1, alternative parameter values
- Explores time-unit alignment, acceleration discretization, temporal decay factor
- Demonstrates A/B testing methodology and parameter sensitivity analysis applied to game AI

This variant is the author's own iteration, kept alongside hzy3.py for comparison of how parameter choices affect tournament performance.

### Layer 3: Strategy Orchestration (`strategy.py`, 546 lines)

Author: Zhenyu He (designed and implemented the high-level strategy orchestration layer).

- Integrates reactive (threat avoidance) + deliberative (goal planning) modes
- Density-adaptive search radius (crowded environments: r=100, sparse: r=800)
- Cost-benefit filtering (do not chase prey smaller than ejection-cost ROI)
- Stateful execution with dynamic abort (collision triggers plan reset)

Architectural inspiration draws from Brooks (1986) subsumption architecture and Boyd OODA loop applied to game AI.

### Layer 4: Game Engine (provided starter)

- `world.py`: game loop, cell initialization, boundary conditions
- `kernel.py`: physics engine (eject, absorb, collision resolution)
- `cell.py`: Cell entity model (position, velocity, radius, ID, collision group)
- `gui.py`: Tkinter visualization
- `consts.py`: tunable constants (world size, frame delta, mass ratio, etc.)

Zhenyu's contribution at this layer: complete reverse-engineering of `kernel.py` physics into the `hzy3.py` rollout simulator (`myeject`, `myabsorb`, `myupdate` helper functions).

### Layer 5: Opponent Library (`sample/` directory, 17 AIs)

- `brownian_motion.py`: weakest baseline (random strategy)
- `cxk.py`: teammate-authored AI with humor elements
- `dynamic.py`: joke AI (sing, dance, rap, basketball)
- `player_ai szx.py`, `player_szx 1.0.py`, `player_szx 1.1.py`: iterative improvements (teammate szx)
- `player.py`: integrated tournament version

Purpose: benchmark gradient for evaluating strategy strength. The AI must beat progressively stronger opponents to be tournament-ready. This layer was co-developed with teammates szx, bs, and others.

## Key Algorithms and Techniques

| Technique | Layer | Purpose |
|---|---|---|
| Relative coordinates + angular enumeration | 1 | Physics-informed action space reduction |
| Torus periodic boundary (9-region expansion) | 1 | Topology-aware search |
| Forward simulation (200-frame rollout) | 1 | Viability verification without heuristics |
| Multi-body collision resolution | 1 (reverse-eng. from Layer 4) | Physics accuracy |
| Reactive + deliberative mix | 3 | Handling unexpected threats gracefully |
| Density-adaptive parameters | 3 | Context-sensitive strategy adjustment |
| Parameter A/B testing | 2 | Scientific tuning methodology |
| Benchmark opponent ladder | 5 | Progressive testing |

## Distinctive Capabilities Demonstrated

1. **Forward simulation for decision verification**: Not pure heuristics; actual physics simulation validates plans before committing.
2. **Physics intuition in algorithm design**: Relative coordinates, torus topology, momentum conservation all embedded in search structure.
3. **Multi-layer architecture**: Reactive safeguards + deliberative planning + orchestration layer.
4. **Reverse-engineering**: Complete mastery of external game engine API (`kernel.py`).
5. **Stateful control**: Sophisticated persistent state machine (not memoryless reaction).
6. **Team collaboration with independent verification**: All team members participated fully in understanding methods before dividing implementation labor.

## Technologies Used

- Python 3
- Object-oriented design (Player, Node, Cell classes)
- Physics simulation
- Game AI architecture (subsumption, OODA)
- Tkinter GUI
- Data structures (linked lists, dictionaries, sets)

## Course Context

- **Institution**: Peking University, School of Earth and Space Sciences
- **Course**: 数据结构与算法 (Python), Data Structures and Algorithms
- **Instructor**: Chen Bin (陈斌)
- **Semester**: Spring 2019 (freshman year)
- **Assignment Type**: Final course project + tournament

## Repository Structure

```
osmo/
├── README.md (this file)
├── CONTRIBUTORS.md
└── -osmo-source0526/
    ├── hzy3.py (Layer 1: primary AI search)
    ├── mechanic.py (Layer 2: parameter variant)
    ├── strategy.py (Layer 3: strategy orchestration)
    └── -osmo-source/
        └── 对战平台代码/ (game engine + opponent library)
            ├── world.py
            ├── kernel.py
            ├── cell.py
            ├── gui.py
            ├── consts.py
            ├── database.py
            └── sample/ (opponent AI library, 17 variants)
```

## License

GNU General Public License v3 (see file headers in `hzy3.py`, `mechanic.py`, `strategy.py`).

## Connections to Later Work

The techniques pioneered in this project informed Zhenyu's subsequent research:

- Forward simulation thinking carried into physics-based constraint satisfaction
- Multi-layer architecture became a template for agent orchestration design
- Parameter sensitivity discipline carried into scientific approach to ML hyperparameter tuning
