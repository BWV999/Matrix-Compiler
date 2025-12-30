# Matrix Calculation Compiler

A domain-specific language (DSL) compiler for matrix operations, built with Flex and Bison. This compiler translates high-level matrix operations into optimized C code.

## Features

### Matrix Operations
- **Basic Arithmetic**: Addition, subtraction, and multiplication
- **Matrix Functions**:
  - `transpose()` - Matrix transposition
  - `det()` - Determinant calculation (supports NxN matrices)
  - `trace()` - Trace calculation
  - `inverse()` - Matrix inversion (supports NxN matrices via Gauss-Jordan elimination)
  - `eigenval()` - Dominant eigenvalue (power iteration method)
  - `eigenvec()` - Dominant eigenvector (power iteration method)

### Advanced Features
- **Arbitrary Matrix Sizes**: Supports NxN square matrices and rectangular matrices
- **Complex Expressions**: Nested operations with parentheses support
- **Automatic Code Generation**: Generates optimized C code with loop unrolling
- **Dimension Checking**: Compile-time validation of matrix dimensions
- **LU Decomposition**: Used internally for determinant calculation

## Requirements

- **Flex** (lexical analyzer generator)
- **Bison** (parser generator)
- **GCC** (C compiler)
- **Math library** (-lm flag)

## Quick Start

### Build and Run

```bash
flex matrix.l
bison -dy matrix.y
gcc lex.yy.c y.tab.c -o matrix_compiler
matrix_compiler < test1.txt
```

Or use the provided script:
```bash
chmod +x build.sh
./build.sh test1.txt
```

## File Structure

```
.
├── matrix.l          # Flex lexer specification
├── matrix.y          # Bison parser and code generator
├── test1.txt         # Comprehensive test suite
├── test2.txt         # Transpose and multiplication tests
├── test3.txt         # Parentheses and complex expressions
├── cmd.txt           # Build commands reference
└── README.md         # This file
```

## Language Syntax

### Matrix Declaration

```c
matrix A(2, 2) = [ 4, 1, 2, 3 ];

matrix B(3, 3) = [ 
  1, 0, 0, 
  0, 1, 0, 
  0, 0, 1 
];
```

### Matrix Operations

```c
// Addition and Subtraction
Sum = A + B;
Diff = A - B;

// Multiplication
Result = A * B;

// Transpose
AT = transpose(A);

// Determinant
det(A);

// Trace
trace(A);

// Inverse
InvA = inverse(A);

// Eigenvalues and Eigenvectors
eigenval(A);
eigenvec(A);
```

### Complex Expressions

```c
// Nested operations
Result = (A + B) * transpose(C);

// Multiple operations
Result = inverse((A + B)) * transpose(C - D);

// Verification
Check = A * inverse(A);  // Should yield identity matrix
```

## Example Usage

### Basic Example (test2.txt)

```c
// 2x3 matrix
matrix A(2, 3) = [
  1, 2, 3,
  4, 5, 6
];

// 3x2 matrix
matrix B(3, 2) = [
  7, 8,
  9, 1,
  2, 3
];

// Multiply: (2x3) * (3x2) = (2x2)
C = A * B;

// Transpose
AT = transpose(A);  // 3x2 matrix
BT = transpose(B);  // 2x3 matrix

// Chain operations: (3x2) * (2x3) = (3x3)
result = AT * BT;
```

### Output

```
A = 
1 2 3 
4 5 6 

B = 
7 8 
9 1 
2 3 

C = 
31 13 
97 43 

transpose(A) = 
1 4 
2 5 
3 6 

transpose(B) = 
7 9 2 
8 1 3 

result = 
31 13 19 
97 43 61 
139 61 87
```

## Advanced Features

### Determinant Calculation (NxN)

The compiler uses LU decomposition for efficient determinant calculation:

```c
matrix E(3, 3) = [
  4, 1, 2,
  1, 5, 3,
  2, 3, 6
];

det(E);  // Calculates determinant
```

### Matrix Inversion (NxN)

Uses Gauss-Jordan elimination with pivoting:

```c
matrix A(4, 4) = [
  5, 2, 1, 1,
  2, 6, 2, 1,
  1, 2, 7, 2,
  1, 1, 2, 8
];

invA = inverse(A);
check = A * invA;  // Verifies identity
```

### Eigenvalue Computation

Power iteration method for dominant eigenvalue/eigenvector:

```c
eigenval(A);  // Computes dominant eigenvalue
eigenvec(A);  // Computes corresponding eigenvector
```

## Implementation Details

### Code Generation

The compiler generates optimized C code with:
- Static array allocation for performance
- Inline loop operations
- Automatic memory management
- Helper functions for complex operations (LU decomposition, Gauss-Jordan, power iteration)

### Error Handling

- Dimension mismatch detection at compile-time
- Singular matrix detection for inversions
- Proper error messages with context

### Limitations

- Maximum matrix size: 100x100 (configurable via MAX_N)
- Eigenvalue/eigenvector: Only computes dominant eigenvalue using power iteration
- Numerical precision: Standard double precision

## Contributing

Contributions are welcome! Areas for improvement:
- Additional matrix operations (QR decomposition, SVD, etc.)
- Sparse matrix support
- Multiple eigenvalue computation
- Matrix factorizations
- Performance optimizations

## License

This project is open source and available under the MIT License.

## Author

Matrix Calculation Compiler - A DSL for Matrix Operations

---

**Note**: This compiler is designed for educational purposes and small to medium-sized matrix computations. For production use with large matrices, consider using specialized libraries like BLAS, LAPACK, or Eigen.
