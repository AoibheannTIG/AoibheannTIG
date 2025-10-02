#!/bin/bash

# Check if circuit file is provided
if [ $# -ne 1 ]; then
    echo "Usage: ./build_circuit.sh <circuit_file.circom>"
    echo "Example: ./build_circuit.sh ../src/magic_square_3x3.circom"
    exit 1
fi

CIRCUIT_FILE=$1
CIRCUIT_NAME=$(basename "$CIRCUIT_FILE" .circom)

# Check if circuit file exists
if [ ! -f "$CIRCUIT_FILE" ]; then
    echo "Error: Circuit file '$CIRCUIT_FILE' not found."
    exit 1
fi

# Create build directory if it doesn't exist
mkdir -p ../build/${CIRCUIT_NAME}

echo "=== Step 1: Compiling the circuit ==="
../../../../../circom/target/release/circom "$CIRCUIT_FILE" --r1cs --wasm --sym --O2 -o ../build/${CIRCUIT_NAME}
if [ $? -ne 0 ]; then 
    echo "Error: Circuit compilation failed"
    exit 1
fi

echo "=== Step 2: Viewing circuit information ==="

# Export R1CS to JSON
snarkjs r1cs export json "../build/${CIRCUIT_NAME}/${CIRCUIT_NAME}.r1cs" "../build/${CIRCUIT_NAME}/${CIRCUIT_NAME}.json"
if [ $? -ne 0 ]; then exit 1; fi



echo -e "\n=== Step 3: Calculating witness ==="
snarkjs wc ../build/${CIRCUIT_NAME}/${CIRCUIT_NAME}_js/${CIRCUIT_NAME}.wasm ../inputs/input_3x3.json ../build/${CIRCUIT_NAME}/${CIRCUIT_NAME}.wtns
if [ $? -ne 0 ]; then 
    echo "Error: Witness calculation failed"
    exit 1
fi
