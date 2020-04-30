======
Theory
======

QR decomposition
----------------

Square matrix
~~~~~~~~~~~~~

Let us consider a square matrix :math:`\mathbf{A}` of size :math:`n`, the QR decomposition of this matrix is:

.. math::
   \mathbf{A} = \mathbf{Q}\mathbf{R}
   :label: QR_decomposition_square

where:

 - :math:`\mathbf{Q}` is an orthogonal square matrix of size :math:`n` such as :math:`\mathbf{Q}^T\mathbf{Q} = \mathbf{Q}\mathbf{Q}^T = \mathbf{I}_n` with :math:`\mathbf{I}_n` the identity matrix of size :math:`n`;
 - :math:`\mathbf{R}` denotes an upper triangular square matrix of size :math:`n`.

The methods ``GetQRDecomposition`` and ``GetFullQRDecomposition`` return the matrices :math:`\mathbf{Q}` and :math:`\mathbf{R}`.

Rectangular matrix
~~~~~~~~~~~~~~~~~~

Let us consider a rectangular matrix :math:`A` of size :math:`m\times n`, the **full** decomposition of this matrix is:

.. math::
   \mathbf{A} = \mathbf{Q}\mathbf{R} = \mathbf{Q}\begin{pmatrix}\mathbf{R}_1 \\ \mathbf{0} \end{pmatrix}
   :label: QR_decomposition_rectangular_full

where:

 - :math:`\mathbf{Q}` is an orthogonal **square** matrix of size :math:`m` such as :math:`\mathbf{Q}^T\mathbf{Q} = \mathbf{Q}\mathbf{Q}^T = \mathbf{I}_n` with :math:`\mathbf{I}_m` the identity matrix of size :math:`n`;
 - :math:`\mathbf{R}` denotes an upper triangular **rectangular** matrix of size :math:`m\times n`;
 - :math:`\mathbf{R}_1` denotes an upper **square** triangular matrix of size :math:`n`;
 - :math:`\mathbf{0}` represents the null **rectangular** matrix of size :math:`(m-n)\times n`.

The matrix :math:`\mathbf{Q}` may be decomposed in two rectangular block matrices :math:`\mathbf{Q}_1` of size :math:`m\times n` and :math:`\mathbf{Q}_2` of size :math:`m\times (m - n)` such as:

.. math::
  \mathbf{Q} = \left(\mathbf{Q}_1, \mathbf{Q}_2\right)

The **thin** decomposition of :math:`\mathbf{A}` is:

.. math::
   \mathbf{A} = \mathbf{Q}_1\mathbf{R}_1
   :label: QR_decomposition_rectangular_thin

For a square matrix, these two decompositions are identical.

In **MathUtils**, the method ``GetQRDecomposition`` applies the thin decomposition and returns the matrices :math:`\mathbf{Q}_1` and :math:`\mathbf{R}_1` whereas the method ``GetFullQRDecomposition`` applies a full decomposition and returns the matrices :math:`\mathbf{Q}` and :math:`\mathbf{R}`.

.. _LS_problem:

Least-square problem
--------------------

Let us consider an overdetermined linear system :math:`\mathbf{A}\mathbf{x} = \mathbf{b}` with :math:`\mathbf{A}` a rectangular matrix of size :math:`m\times n` and :math:`\mathbf{b}` a vector of size :math:`m`. We assume the following condition:

.. math::
   m > n

The solving of this system is achieved in **MathUtils** by using a least-square method (``LeastSquareSolver``), based on the **Eigen** method: ``bdcSvd``:

.. math::
   \min\left\Vert\mathbf{A}\mathbf{x} - \mathbf{b}\right\Vert^2

The solution :math:`\mathbf{x}` satisfies:

.. math::
   \mathbf{x} = \mathbf{A}^{-1}\mathbf{b}
   :label: Solution_LS_problem

where :math:`\mathbf{A}^{-1}` is the pseudoinverse matrix of :math:`\mathbf{A}`.

Least-square problem subject to equality constraints
----------------------------------------------------

Let us consider the same problem as presented in :ref:`LS_problem` but we now add an equality constraint, so the problem becomes:

.. math::
   \min\left\Vert\mathbf{A}\mathbf{x} - \mathbf{b}\right\Vert^2 \: \text{subject to} \: \mathbf{C}\mathbf{x} = \mathbf{d}

with :math:`\mathbf{C}` a rectangular matrix of size :math:`p\times n` and :math:`\mathbf{d}` a vector of size :math:`p`.

For solving this problem, we follow the demonstration of [Grivet2016]_ in the section D.3.2. We start by doing a QR decomposition to the matrix :math:`\mathbf{C}^T` using :eq:`QR_decomposition_rectangular_full`:

.. math::
   \mathbf{C}^T = \mathbf{Q}\begin{pmatrix}\mathbf{R} \\ \mathbf{0} \end{pmatrix}

with :math:`\mathbf{Q}` a square matrix of size :math:`n` and :math:`\mathbf{R}` a square matrix of size :math:`p`.

Let us define the vector :math:`\mathbf{y}` of size :math:`p` and the vector :math:`\mathbf{z}` of size :math:`n-p` such as:

.. math::
   \mathbf{x} = \mathbf{Q}\begin{pmatrix} \mathbf{y} \\ \mathbf{z} \end{pmatrix}
   :label: Definition_y_z

By using the orthogonality of the matrix :math:`\mathbf{Q}`, we get:

.. math::
   \mathbf{C}\mathbf{x} = \mathbf{R}^T\mathbf{y}

Thus, the unknown :math:`\mathbf{y}` is determined from:

.. math::
   \mathbf{R}^T\mathbf{y} = \mathbf{d}
   :label: Linear_system_y

Let us define the two rectangular block matrices :math:`\mathbf{A}_1` of size :math:`m\times p` and :math:`\mathbf{A}_2` of size :math:`m\times (n - p)` such as:

.. math::
  \mathbf{A}\mathbf{Q} = \left(\mathbf{A}_1, \mathbf{A}_2\right)

We have:

.. math::
   \mathbf{A}\mathbf{x} = \mathbf{A}_1\mathbf{y} + \mathbf{A}_2\mathbf{z}

Using :eq:`Linear_system_y`, the unknown :math:`\mathbf{z}` is figured out by the least-square problem:

.. math::
   \min\left\Vert\mathbf{A}_2\mathbf{z} - \left(\mathbf{b} - \mathbf{A}_1\mathbf{R}^{-T}\mathbf{d}\right)\right\Vert^2

Finally, the solution :math:`\mathbf{x}` is obtained from :eq:`Definition_y_z`.

This solution satisfies exactly the equality constraint :math:`\mathbf{C}\mathbf{x} = \mathbf{d}` but not exactly the linear system :math:`\mathbf{A}\mathbf{x} = \mathbf{b}`. Thus, contrairy to :ref:`LS_problem`, the solution does not satisfy :eq:`Solution_LS_problem`. Some errors appear because of the equality constraint [Grivet2016]_.

This method is used in ``LeastSquareSolverConstraint``.

.. [Grivet2016] S. Grivet-Talocia and B. Gustavsen. Passive macromodeling. Theory and applications. 2016.