# ZK Playground

A zero-knowledge proof playground containing implementations and optimization tools for zero-knowledge circuits.

## Project Structure

```
ZK_playground/
├── magic_square/                    # Magic Square ZK Proof Implementation
│   ├── src/                        # Source code files
│   │   ├── magic_square_3x3.circom # Circom circuit definition
│   │   └── magic_square_validator.py # Python validator
│   ├── inputs/                     # Input data files
│   │   └── input_3x3.json         # Sample 3x3 magic square input
│   ├── scripts/                    # Build scripts
│   │   └── build_circuit.sh       # Circuit compilation script
│   ├── build/                      # Build artifacts (generated)
│   └── README.md                   # Magic square documentation
├── r1cs_optimize.py               # R1CS optimization tool
└── README.md                      # This file
```

## Components

### Magic Square Implementation
A zero-knowledge proof system for verifying 3x3 magic squares using Circom. The circuit verifies that a 3x3 grid is a magic square by checking:
- All rows, columns, and diagonals sum to the magic constant
- The sum of all numbers equals 45 (1+2+...+9)
- The sum of squares equals 285 (1²+2²+...+9²)

See the [magic_square/README.md](magic_square/README.md) for detailed usage instructions.

### R1CS Optimization Tool
A Python tool (`r1cs_optimize.py`) that optimizes R1CS (Rank-1 Constraint Systems) by:
1. Eliminating internal variables iteratively
2. Applying algebraic factoring to simplify structure
3. Reducing constraints to the form a*b + c = 0

## Usage

### Building the Magic Square Circuit
```bash
cd magic_square/scripts
./build_circuit.sh ../src/magic_square_3x3.circom
```

### Optimizing R1CS
```bash
python r1cs_optimize.py magic_square/build/magic_square_3x3/magic_square_3x3.json -o magic_square/build/magic_square_3x3/magic_square_3x3_optimised.json
```

### Validating Magic Squares
```bash
cd magic_square/src
python magic_square_validator.py
```

## Requirements

- Python 3.x
- Circom compiler
- Node.js (for witness generation)
- SymPy (for R1CS optimization)

## Getting Started

1. Navigate to the magic square directory: `cd magic_square`
2. Follow the instructions in the magic square README to build the circuit
3. Use the R1CS optimization tool to optimize the generated constraints
4. Generate witnesses and test the zero-knowledge proof system