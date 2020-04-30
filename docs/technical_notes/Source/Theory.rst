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

The matrix :math:`\mathbf{Q}` may be decomposed in two block matrices :math:`\mathbf{Q}_1` and :math:`\mathbf{Q}_2` such as:

.. math::
  \mathbf{Q} = \left(\mathbf{Q}_1, \mathbf{Q}_2\right)

:math:`\mathbf{Q}_1` and :math:`\mathbf{Q}_2` are two **rectangular** matrices of size :math:`m\times n` and :math:`m\times (m - n)`, respectively.

The **thin** decomposition of :math:`\mathbf{A}` is:

.. math::
   \mathbf{A} = \mathbf{Q}_1\mathbf{R}_1
   :label: QR_decomposition_rectangular_thin

For a square matrix, these two decompositions are identical.

In **MathUtils**, the method ``GetQRDecomposition`` applies the thin decomposition and returns the matrices :math:`\mathbf{Q}_1` and :math:`\mathbf{R}_1` whereas the method ``GetFullQRDecomposition`` applies a full decomposition and returns the matrices :math:`\mathbf{Q}` and :math:`\mathbf{R}`.
