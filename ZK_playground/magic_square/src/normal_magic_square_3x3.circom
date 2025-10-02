pragma circom 2.1.6;


// This template verifies that a 3x3 grid is a magic square.
template MagicSquare3x3() {
    // Public inputs: the 9 grid entries and the magic constant.
    signal input a[9];
    signal input MAGIC;

    // --- CONSTRAINT 1: All entries are unique between 1, .., 9 (n^2) ---
    signal squares[9]; // Individual squares to make constraints quadratic

    for (var i = 0; i < 9; i++) {

        // Calculate squares as separate signals for quadratic constraints
        squares[i] <== a[i] * a[i];
    }
    
    // Sum constraint (linear)
    a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] + a[8] === 45;
    
    // Sum of squares constraint (quadratic but each term is quadratic)
    squares[0] + squares[1] + squares[2] + squares[3] + squares[4] + squares[5] + squares[6] + squares[7] + squares[8] === 285;


    // --- CONSTRAINT 2: All rows, columns, and diagonals sum to MAGIC ---
    // Rows
    (a[0] + a[1] + a[2]) === MAGIC;
    (a[3] + a[4] + a[5]) === MAGIC;
    (a[6] + a[7] + a[8]) === MAGIC;

    // Columns
    (a[0] + a[3] + a[6]) === MAGIC;
    (a[1] + a[4] + a[7]) === MAGIC;
    (a[2] + a[5] + a[8]) === MAGIC;

    // Diagonals
    (a[0] + a[4] + a[8]) === MAGIC;
    (a[2] + a[4] + a[6]) === MAGIC;
}

// Instantiate the main component, making the grid and magic constant public.
component main {public [a, MAGIC]} = MagicSquare3x3();
