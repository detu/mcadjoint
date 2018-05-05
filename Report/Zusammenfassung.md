Im Moment arbeite ich daran, dass der Monte-Carlo-Schätzer konvergiert.
Dazu muss die Neumann-Reihe einer Matrix

$$ A := \begin{pmatrix} I - J_{00}^T  & -J_{10}^T & & & \\  & I-J_{11}^T &  -J_{21}^T & & & & \\ \text{usw}\end{pmatrix} $$

konvergieren. Dabei sind

$$ J_{ii} := \begin{pmatrix} \frac{\partial \Pi^{(i)}}{\partial p^{(i)}} &0 \\\frac{\partial \Sigma^{(i)}}{\partial p^{(i)}} & I\end{pmatrix}$$

und $$J_{i, i-1} := \begin{pmatrix} 0 & \frac{\partial \Pi^{(i)}}{\partial S^{(i-1)}}{}\\ 0 & \frac{\partial\Sigma^{(i)}}{\partial S^{(i-1)}}\end{pmatrix} $$

Der hochgestellte Index $(i)$ steht für einen Zeitschritt.
$\Pi$ ist der Vektor der Residuen der Druckgleichung und $\Sigma$ der Vektor der Residuen der Sättigungsgleichung. $S$ ist der Vektor der Sättigungen und $p$ der Vektor der Drücke (Siehe report.pdf).

Aus $A$ kann man eine Übergangsmatrix für ein Markovkettenverfahren gewinnen. Laut dem Paper ist es in diesem Fall ziemlich gut, einfach A zeilenstochastisch zu machen.

Das habe ich so gemacht. Das Problem ist nun, dass $\frac{\partial \Pi^{(i)}}{\partial p^{(i)}}$ sehr grosse Einträge hat im Vergleich zu den anderen Blöcken. Die Random Walks über die Indices bleiben also in den Druckzuständen stecken und kommen nicht weiter.

Das hat zwei Konsequenzen:
- Da $\frac{\partial \Pi^{(i)}}{\partial p^{(i)}}$ auf der Blockdiagonale der Matrix $A$ steht, kann ich das Vorwärtsproblem nicht einen Schritt weiterlaufen lassen, sonder muss warten, bis durch Zufall ein Index 
