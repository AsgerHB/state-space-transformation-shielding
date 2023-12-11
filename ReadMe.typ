#set math.equation(numbering: "(1)")

= Read Me

This set of experiments was created from the idea that observations computed from physical equations might provide a superior signal to learning algorithms, compared to the raw state of the system. For example, the total mechanical energy contained in a physical system is

$ E_"mek" = E_"pot" + E_"kin" = m g h + 1/2 m v^2 $ <EMek>

which might work better for learning a strategy for the bouncing ball, compared to just observing $v$ and $p$.

== Experiments on the Bouncing Ball - story

Training a strategy to hit the bouncing ball by observing the full state $(v, p)$ yields an average performance of 40 swings per 120 seconds. Unfortunately, observing only $E_"mek"$ uses 146 swings per 120 seconds. This is because the outcome of hitting the ball depends on $v$ and $p$, which leads to a lot of useless swings from the strategy. 

Observing $E_"mek"$ in addition to $(v, p)$ leads to [data missing]

Since the hit#footnote([This is slightly confusing because choosing the "hit" action does not necesarily hit the ball. It just swings at it. Only if $-4 < v and p > 4$ does the ball actually get hit. See @PossibleToHit]) will only connect for certain values of $v$ and $p$, I chose to also count the number of times the hits actually connected. This meant that the strategy observing just $E_"mek"$ would not be penalised for all the unnecesary swings that had no effect. By this criterion, $E_"mek"$ strategy then achieved a decent performance of 41 hits per 120 seconds. Which is close to the strategy observing $(v, p)$. 

Obseriving both $(v, p)$ and $E_"mek"$ led to 37 swings per 120 seconds. 

Running the queries once, there does not seem to be any significant difference between the learning outcomes, except that $E_"mek"$ alone is not enough to efficiently hit the ball, if useless swings are counted.

The queries would have to be repeated multiple times with a set number of training runs to discover if there is any difference.q

== Experiments on the Bouncing Ball

The full set of learning queries is seen in @BBQueries. 

Strategies were learned using observations of either the full state, $(v, p)$, or the mechanical energy formula in @EMek. Or some combination of the two.

The learner was either trained to minimize the number of "swings" or the number of "hits."  This is because the mechanical energy is not enough to determine whether it is possible to hit the ball at the given state. This forces the strategy to take a lot of unnecessary swings in order to be sure to hit the ball. Looking instead only at how often the ball was hit may be more "fair."

The result of running the experiment once with default learning parameters shows 40 swings and 37 hits, when trying to minimize the number of swings. Minimizing just hits, leads to 497 swings and 39 hits.

Observing only $E_"mek"$ Leads to 146 swings and 44 hits when trying to minimize swings. Minimizing for hits seems to be the same, with 41 hits.

Also observing $p$ gets the value down to 41 or 38, depending on if it is minimizing for swings or for hits.

== Visualisations for the Bouncing Ball

Visualisations showing the mechanical energy of the bouncing ball for different states $(v, p)$ are shown in @MechanicalEnergySmooth and @MechanicalEnergyGrouped. 

#table(columns:2, stroke: none, 
  [#figure(image("BB Mechanical Energy Smooth.png"),
    caption:[Bouncing ball $E_"mek"$ values.],
    ) <MechanicalEnergySmooth>],

  [#figure(image("BB Mechanical Energy Grouped.png"),
    caption:[Bouncing ball $E_"mek"$ values grouped into 10 discrete categories.]
    ) <MechanicalEnergyGrouped>],
  [#figure(image("BBShield Overlaid Mechanical Energy.png"),
  caption:[Bouncing ball shield overlaid with @MechanicalEnergyGrouped.]
  ) <BBShieldOverlaidMechanicalEnergy>],
  [#figure(image("Possible to Hit.png"),
  caption:[Overlay to @MechanicalEnergyGrouped showing states where it is possible to hit the ball.]
  ) <PossibleToHit>]
)

It is not very scientific, but I made sure the axis matched a visualisation of the safety strategy for the bouncing ball, and then I used a graphics editing program to overlay one with the other, which is shown in @BBShieldOverlaidMechanicalEnergy. Thus, the visualisation is not guaranteed to be pixel perfect, but it is very nearly correct.

#pagebreak()

= Appendix

#figure(image("BB Queries.png"), caption:[Learning queries.])<BBQueries>