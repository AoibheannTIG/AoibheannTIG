# Magic Square 3x3 Zero-Knowledge Proof

This project implements a zero-knowledge proof system for verifying 3x3 magic squares using Circom.

## Project Structure

```
magic_square/
├── src/                          # Source code files
│   ├── magic_square_3x3.circom   # Circom circuit definition
│   └── magic_square_validator.py # Python validator for magic squares
├── inputs/                       # Input data files
│   └── input_3x3.json           # Sample input for 3x3 magic square
├── scripts/                      # Build and utility scripts
│   └── build_circuit.sh          # Script to compile and build the circuit
├── build/                        # Build artifacts (generated)
└── README.md                     # This file
```

## Files Description

### Source Files (`src/`)
- **`magic_square_3x3.circom`**: The main Circom circuit that defines the constraints for verifying a 3x3 magic square
- **`magic_square_validator.py`**: Python utility to validate magic squares independently

### Input Files (`inputs/`)
- **`input_3x3.json`**: Sample input data containing a valid 3x3 magic square and its magic constant

### Scripts (`scripts/`)
- **`build_circuit.sh`**: Build script that compiles the Circom circuit and generates witness

## Usage

### Building the Circuit

```bash
cd scripts
./build_circuit.sh ../src/magic_square_3x3.circom
```

### Validating Magic Squares (Python)

```bash
cd src
python magic_square_validator.py
```

## Circuit Details

The Circom circuit verifies that a 3x3 grid is a magic square by checking:
1. All rows sum to the magic constant
2. All columns sum to the magic constant  
3. Both diagonals sum to the magic constant
4. The sum of all numbers equals 45 (1+2+...+9)
5. The sum of squares equals 285 (1²+2²+...+9²)

## Input Format

The input JSON file should contain:
```json
{
  "a": [2, 7, 6, 9, 5, 1, 4, 3, 8],
  "MAGIC": 15
}
```

Where `a` is the flattened 3x3 grid (row by row) and `MAGIC` is the magic constant.
