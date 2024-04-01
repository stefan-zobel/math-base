/*
 * Copyright 2024 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package math.linalg;

import java.util.Arrays;

import math.gemm.Dgemm;
import math.gemm.Trans;

/**
 * A minimal matrix of doubles that can't do much more than addition,
 * subtraction and multiplication.
 */
public class DMatrix {

    protected final int rows;
    protected final int cols;
    protected final double[] a;

    public static DMatrix identity(int dim) {
        DMatrix I = new DMatrix(dim, dim);
        for (int i = 0; i < dim; ++i) {
            I.setUnsafe(i, i, 1.0);
        }
        return I;
    }

    public DMatrix(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        this.a = new double[checkArrayLength(rows, cols)];
    }

    public DMatrix(DMatrix other) {
        this(other.rows, other.cols, Arrays.copyOf(other.a, other.a.length));
    }

    protected DMatrix(int rows, int cols, double[] a) {
        this.rows = rows;
        this.cols = cols;
        this.a = a;
    }

    public DMatrix copy() {
        return new DMatrix(this);
    }

    public int numColumns() {
        return cols;
    }

    public int numRows() {
        return rows;
    }

    public double get(int row, int col) {
        checkIndex(row, col);
        return a[idx(row, col)];
    }

    public double getUnsafe(int row, int col) {
        return a[idx(row, col)];
    }

    public DMatrix set(int row, int col, double val) {
        checkIndex(row, col);
        a[idx(row, col)] = val;
        return this;
    }

    public void setUnsafe(int row, int col, double val) {
        a[idx(row, col)] = val;
    }

    public double[] getArrayUnsafe() {
        return a;
    }

    public DMatrix scale(double alpha) {
        DMatrix B = new DMatrix(rows, cols);
        if (alpha == 0.0) {
            return B;
        }
        double[] _a = a;
        double[] _b = B.a;
        for (int i = 0; i < _b.length; ++i) {
            _b[i] = alpha * _a[i];
        }
        return B;
    }

    public DMatrix transpose() {
        if (rows == 1 || cols == 1) {
            return new DMatrix(cols, rows, Arrays.copyOf(a, a.length));
        }
        DMatrix AT = new DMatrix(cols, rows);
        int cols_ = cols;
        int rows_ = rows;
        for (int col = 0; col < cols_; ++col) {
            for (int row = 0; row < rows_; ++row) {
                AT.setUnsafe(col, row, getUnsafe(row, col));
            }
        }
        return AT;
    }

    public DMatrix add(DMatrix B) {
        checkEqualDimension(this, B);
        DMatrix C = new DMatrix(rows, cols);
        double[] _a = a;
        double[] _b = B.a;
        double[] _c = C.a;
        for (int i = 0; i < _a.length; ++i) {
            _c[i] = _a[i] + _b[i];
        }
        return C;
    }

    public DMatrix minus(DMatrix B) {
        checkEqualDimension(this, B);
        DMatrix C = new DMatrix(rows, cols);
        double[] _a = a;
        double[] _b = B.a;
        double[] _c = C.a;
        for (int i = 0; i < _a.length; ++i) {
            _c[i] = _a[i] - _b[i];
        }
        return C;
    }

    public DMatrix mul(DMatrix B) {
        checkMul(this, B);
        DMatrix C = new DMatrix(this.rows, B.cols);
        Dgemm.dgemm(Trans.NO_TRANS, Trans.NO_TRANS, C.rows, C.cols, cols, 1.0, a, 0, rows, B.a, 0, B.rows, 0.0, C.a, 0,
                C.rows);
        return C;
    }

    protected final int idx(int row, int col) {
        return col * rows + row;
    }

    protected void checkIndex(int row, int col) {
        if (row < 0 || row >= rows) {
            throw new IllegalArgumentException("Illegal row index " + row + " in (" + rows + " x " + cols + ") matrix");
        }
        if (col < 0 || col >= cols) {
            throw new IllegalArgumentException(
                    "Illegal column index " + col + " in (" + rows + " x " + cols + ") matrix");
        }
    }

    protected static void checkSameRows(DMatrix A, DMatrix B) {
        if (A.numRows() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != B.numRows() (" + A.numRows() + " != " + B.numRows() + ")");
        }
    }

    protected static void checkSameCols(DMatrix A, DMatrix B) {
        if (A.numColumns() != B.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numColumns() (" + A.numColumns() + " != " + B.numColumns() + ")");
        }
    }

    protected static void checkEqualDimension(DMatrix A, DMatrix B) {
        checkSameRows(A, B);
        checkSameCols(A, B);
    }

    protected static void checkMul(DMatrix A, DMatrix B) {
        if (A.numColumns() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numRows() (" + A.numColumns() + " != " + B.numRows() + ")");
        }
    }

    protected static void checkMul(DMatrix A, DMatrix B, DMatrix C) {
        if (A.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numRows() != C.numRows() (" + A.numRows() + " != " + C.numRows() + ")");
        }
        if (A.numColumns() != B.numRows()) {
            throw new IndexOutOfBoundsException(
                    "A.numColumns() != B.numRows() (" + A.numColumns() + " != " + B.numRows() + ")");
        }
        if (B.numColumns() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numColumns() != C.numColumns() (" + B.numColumns() + " != " + C.numColumns() + ")");
        }
    }

    protected static void checkAdd(DMatrix A, DMatrix B, DMatrix C) {
        checkEqualDimension(A, B);
        if (B.numRows() != C.numRows()) {
            throw new IndexOutOfBoundsException(
                    "B.numRows() != C.numRows() (" + B.numRows() + " != " + C.numRows() + ")");
        }
        if (B.numColumns() != C.numColumns()) {
            throw new IndexOutOfBoundsException(
                    "B.numColumns() != C.numColumns() (" + B.numColumns() + " != " + C.numColumns() + ")");
        }
    }

    protected static int checkArrayLength(int rows, int cols) {
        long length = (long) checkRows(rows) * (long) checkCols(cols);
        if (length > (long) Integer.MAX_VALUE) {
            throw new IllegalArgumentException(
                    "rows x cols (= " + length + ") exceeds the maximal possible length (= 2147483647) of an array");
        }
        return (int) length;
    }

    protected static int checkRows(int rows) {
        if (rows <= 0) {
            throw new IllegalArgumentException("number of rows must be strictly positive : " + rows);
        }
        return rows;
    }

    protected static int checkCols(int cols) {
        if (cols <= 0) {
            throw new IllegalArgumentException("number of columns must be strictly positive : " + cols);
        }
        return cols;
    }
}
