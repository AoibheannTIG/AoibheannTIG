def is_magic_square(square):
    """
    Verifies if a 2D list is a magic square.

    A magic square is a square grid of unique numbers between 1, .., n^2 where the sum of each row, column,
    and both main diagonals is the same. 

    Args:
        square: A 2D list of numbers representing the square.

    Returns:
        True if it is a magic square, False otherwise.
    """
    # First, check if the input is actually a square (n x n).
    n = len(square)
    for row in square:
        if len(row) != n:
            print("This is not a square matrix.")
            return False

    # Check if every input is a unique integer between 1,.., n^2.

    sums = 0
    sqs = 0

    for row in square:
        for num in row:
            sums += num
            sqs += num**2
    n_sq = n**2
    target_sum = n_sq*(n_sq + 1)/2
    target_sqs = n_sq*(n_sq + 1)*(2*n_sq + 1)/6
    if sums != target_sum:
        print(f"The sum of the inputs {sums} must be {target_sum}.")
        return False
    if sqs != target_sqs:
        print(f"The sum of the squares of the inputs {sqs} must be {target_sqs}.")
        return False

    # Calculate the magic constant from the first row's sum.
    magic_constant = sum(square[0])
    print(f"The expected magic constant is: {magic_constant}")

    # 1. Check the sums of all rows.
    print("\nChecking rows...")
    for i in range(n):
        row_sum = sum(square[i])
        print(f"Sum of row {i+1}: {row_sum}")
        if row_sum != magic_constant:
            print(f"Row {i+1} sum does not match the magic constant.")
            return False

    # 2. Check the sums of all columns.
    print("\nChecking columns...")
    for j in range(n):
        col_sum = 0
        for i in range(n):
            col_sum += square[i][j]
        print(f"Sum of column {j+1}: {col_sum}")
        if col_sum != magic_constant:
            print(f"Column {j+1} sum does not match the magic constant.")
            return False

    # 3. Check the sum of the main diagonal (top-left to bottom-right).
    print("\nChecking main diagonal...")
    main_diag_sum = 0
    for i in range(n):
        main_diag_sum += square[i][i]
    print(f"Sum of main diagonal: {main_diag_sum}")
    if main_diag_sum != magic_constant:
        print("Main diagonal sum does not match the magic constant.")
        return False

    # 4. Check the sum of the anti-diagonal (top-right to bottom-left).
    print("\nChecking anti-diagonal...")
    anti_diag_sum = 0
    for i in range(n):
        anti_diag_sum += square[i][n - 1 - i]
    print(f"Sum of anti-diagonal: {anti_diag_sum}")
    if anti_diag_sum != magic_constant:
        print("Anti-diagonal sum does not match the magic constant.")
        return False

    # If all checks passed, it's a magic square!
    return True

# --- Example Usage ---

# Example 1: A valid 3x3 magic square
print("--- Verifying Example 1 (A valid magic square) ---")
magic_square_example = [
    [2, 7, 6],
    [9, 5, 1],
    [4, 3, 8]
]

if is_magic_square(magic_square_example):
    print("\nResult: Congratulations! This is a magic square.")
else:
    print("\nResult: This is not a magic square.")

print("\n" + "="*50 + "\n")

# Example 2: An invalid 4x4 square
print("--- Verifying Example 2 (An invalid square) ---")
not_magic_square_example = [
    [16, 3, 2, 13],
    [5, 10, 11, 8],
    [9, 6, 7, 12],
    [4, 15, 14, 2]  # Changed the last element from 1 to 2 to make it invalid
]

if is_magic_square(not_magic_square_example):
    print("\nResult: Congratulations! This is a magic square.")
else:
    print("\nResult: This is not a magic square.")
