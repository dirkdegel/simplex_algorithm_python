# Simplex Algorithm Python Code
Simplex algorithm in python code for linear programming. 

This implementation is inteded to show the general ideas of the simplex algorithm as proposed by **George Dantzig** in 1947. It is **not** suitable for productive use and may include errors.

The algorithm solves linear problems in canonical form:
$$ \max\{ c^\top x \mid Ax \leq b, x \geq 0, \} $$
where $A\in \mathbb{R}^{m \times n}, b \in \mathbb{R}^{m}, b\geq 0$ and $c \in \mathbb{R}^{n}$.

The algorithm first constructs a **tableau** of the following form:
$$
\begin{array}{|c|cc|c|}
\hline
1 & -c^\top & \mathbf{0} & 0 \\ \hline
\mathbf{0} & A & \mathbf{I}_m & b \\ \hline
\end{array}
$$
Then iteratively new tableaus are generated until the optimal solution is found or it is proven that $c^\top x$ is unbounded over $\{Ax \leq b, x \geq 0\}$.
