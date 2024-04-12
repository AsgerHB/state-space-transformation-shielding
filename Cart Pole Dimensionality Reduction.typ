= Cart Pole Dimensionality Reduction
Recall that the state space for the cartpole is 
$(x, dot(x), θ, dot(θ))^top$.

Alrighty so contrary to the bouncing ball, there is no obvious invariant that we could think of for the cart pole setup.
And while it is easy to generate a shield that keeps *either* the cart's position or the pole's angle within bounds, synthesising a shield that does both seems to be infeasible. Yesterday evening I started a 3-million partition grid, which took 1.2 hours to complete, and even then did not yield a safe strategy. Just another black square. 

Safety strategies could be created in decent resolution for just the position or angle, it's just that with the 4 dimensions together, it becomes infeasible. The safety strategy for position alone is shown in @position_regular, and the safety strategy for the angle in @angle_regular. These shields are of size $150 times 150$.

#grid(columns: 2, 
  [#figure(image("Graphics/Cart Pole/Position Regular.png", fit:"cover"), 
      caption: "Shield for position only.") 
    <position_regular>],

  [#figure(image("Graphics/Cart Pole/Angle Regular.png", fit:"cover"), 
      caption: "Shield for angle only.") 
    <angle_regular>],

)

The hope is that we can use knowledge gained from these two independent shields, to obtain an abstraction that shields the whole system.
For each shield, we can see some similarity in the shape of the upper and lower decision boundaries. This is not surprising, since the problem itself is symmetrical. If nothing else, the decision boundary is skewed compared to the alignment of the partitioning. So these are shapes that are not so easy to capture using squares.

The idea then was to fit a polynomial to the decision boundaries. @position_polynomial_fitting_top shows a polynomial fit to the topmost decision boundary. #footnote[I wrote a function that took a 2D-grid and two values, and returned samples from the middle of partitions that had the first values and which bordered a partition that had the second value.] Since the two boundaries aren't exactly alike, I used the average values of the two decision boundaries to fit a polynomial, as shown in @position_polynomial_fitting_average.

#grid(columns: 2, 
  [#figure(image("Graphics/Cart Pole/Position Polynomial Fitting Top.png", fit:"cover"), 
      caption: "Fitting a polynomial to the topmost decision border of the position-shield.") 
    <position_polynomial_fitting_top>],

  [#figure(image("Graphics/Cart Pole/Position Polynomial Fitting Average.png", fit:"cover"), 
      caption: "Fitting a polynomial to the average of the two decision borders") 
    <position_polynomial_fitting_average>],

)

So now I have a polynomial $P$, where $P(x) = dot(x)$ such that $dot(x)$ is at the decision boundary. An offset $omega$ can be applied, such that the boundary can be shifted: $P(x) - omega = dot(x)$ which can be re-written as: $omega = P(x) - dot(x)$. 
Now, creating a grid whose state space is $vec(x, omega)$ instead of $vec(x, dot(x))$ should create a decision boundary that better aligns with the partition boundary.

This is an abstraction for which we can still compute reachability, since $dot(x)$ is given from $x$ and $omega$. Let $P_1(x, dot(x)) = P(x) - dot(x)$. Then $P_1^(-1)(x, omega) = P(x) - omega = dot(x)$. A safety strategy in the $vec(x, omega)$-grid abstraction can be used as a safety strategy to keep the cart's position within bounds, by looking up the set of allowed actions as $vec(x, P_1(x))$. 

@position_10th_degree_polynomial shows a 10th degree polynomial, fit to the average of the two decision boundaries. 
@position_1st_degree_polynomial shows the same for a 1st degree polymonial. Note that neither decision boundary is completely straight, because the upper and lower boundary aren't the same polynomial. However, This is still a better abstraction.
Using this altered state space, I was able to shield the cart's position at much coarser granularity than in the original state space. 
@position_1st_degree_polynomial_coarse and @position_1st_degree_polynomial_very_coarse show the shields at coarser granularities, where a straight decision boundary (though more restrictive) can be successfully used.

The same thing can be seen for the shield that controls the pole's angle in @angle_10th_degree_polynomial through @angle_1st_degree_polynomial_very_coarse. Incidentally, it seems that @angle_1st_degree_polynomial_very_coarse could be a safety strategy for the system as a whole, since it is so restrictive in keeping the pole upright, that the cart hardly seems to move. However, it might be the case that it will eventually allow the cart to drift off to one side.

Unfortunately, combining the two polynomials into a state space $(x, P_1(x, dot(x)), θ, P_2(θ, dot(θ)))^top$ 
does not yield a safety strategy at any granularity that I have tried. 

#grid(columns: 2,
  [#figure(image("Graphics/Cart Pole/Position 10th Degree Polynomial.png", fit:"cover"), 
      caption: "Fitting a polynomial to the average of the two decision borders") 
    <position_10th_degree_polynomial>],

  [#figure(image("Graphics/Cart Pole/Position 1st Degree Polynomial.png", fit:"cover"), 
      caption: "Fitting a polynomial to the topmost decision border of the  position-shield.") 
    <position_1st_degree_polynomial>],

  [#figure(image("Graphics/Cart Pole/Position 1st Degree Polynomial Coarse.png", fit:"cover"), 
      caption: "Fitting a polynomial to the average of the two decision borders") 
    <position_1st_degree_polynomial_coarse>],

  [#figure(image("Graphics/Cart Pole/Position 1st Degree Polynomial Very Coarse.png", fit:"cover"), 
      caption: "Fitting a polynomial to the topmost decision border of the position-shield.") 
    <position_1st_degree_polynomial_very_coarse>],

)

Future work would be to find better expressions for the decision boundaries, perhaps by combining the polynomials depending on the sign of $dot(x)$. 

#grid(columns: 2,
  [#figure(image("Graphics/Cart Pole/Angle 10th Degree Polynomial.png", fit:"cover"), 
      caption: "Fitting a polynomial to the average of the two decision borders") 
    <angle_10th_degree_polynomial>],

  [#figure(image("Graphics/Cart Pole/Angle 1st Degree Polynomial.png", fit:"cover"), 
      caption: "Fitting a polynomial to the topmost decision border of the position-shield.") 
    <angle_1st_degree_polynomial>],

  [#figure(image("Graphics/Cart Pole/Angle 1st Degree Polynomial Coarse.png", fit:"cover"), 
      caption: "Fitting a polynomial to the average of the two decision borders") 
    <angle_1st_degree_polynomial_coarse>],

  [#figure(image("Graphics/Cart Pole/Angle 1st Degree Polynomial Very Coarse.png", fit:"cover"), 
      caption: "Fitting a polynomial to the topmost decision border of the position-shield.") 
    <angle_1st_degree_polynomial_very_coarse>],

)