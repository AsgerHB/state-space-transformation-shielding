#set math.equation(numbering: "(1)")

= Read Me

This set of experiments was created from the idea that observations computed from physical equations might provide a superior signal to learning algorithms, compared to the raw state of the system. For example, the total mechanical energy contained in a physical system is

#let mek = "mek"
#let pot = "pot"
#let kin = "kin"

$ E_mek = E_pot + E_kin = m g h + 1/2 m v^2 $ <EMek>

which might work better for learning a strategy for the bouncing ball, compared to just observing $v$ and $p$.

== Experiments on the Bouncing Ball

Note the difference between swings and hits. The agent can swing at the ball at any time, but will only hit the ball if $4 < v and -4 < p$. This area is illustated in @PossibleToHit. Both values are included here because not all sets of observatiotns make it possible to deduce whether a swing will actually be hit in the given state.

The data from @tab:LearningOutcomes can be found in the appendix, @BBQueries.
The table shows learning outcomes from strategies trained with different sets of observations. Each strategy was trained once with 3000 episodes.

The observation $Delta E_mek$ represents the mechanical energy gained from hitting the ball in the current state. This value is computed by the function `delta_e_mek`.
 
#figure(caption: [Learning outcomes based on different observations. Scores are expected outcome of a 120 second trace.],
  table(columns: 3,
    strong[Observations], strong[Swings], strong[Hits],
    $v, p$, [43], [39],
    $E_mek$, [157], [44],
    $E_mek, p$, [202], [42],
    $E_mek, p, v$, [41], [40],
    $E_mek, Delta E_mek$, [39], [37],
    $E_kin, E_pot$, [214], [40],
  )
)<tab:LearningOutcomes>

NB: Since training is stochastic, a difference from e.g. 38 to 44 is not statistically significant.

Observing just $E_mek$ is not enough to determine if a swing will hit, and so the performance suffers. Observing $E_mek$ in addition to $p$ and $v$ has no significant impact on performance. Observing $E_mek$ and $Delta E_mek$ is almost as good as observing $v$ and $p$ directly.

== Visualisations for the Bouncing Ball

Visualisations showing the mechanical energy of the bouncing ball for different states $(v, p)$ are shown in @MechanicalEnergySmooth and @MechanicalEnergyGrouped. 

#table(columns:2, stroke: none, 
  [#figure(image("Graphics/BB Mechanical Energy Smooth.png"),
    caption:[Bouncing ball $E_mek$ values.],
    ) <MechanicalEnergySmooth>],

  [#figure(image("Graphics/BB Mechanical Energy Grouped.png"),
    caption:[Bouncing ball $E_mek$ values grouped into 10 discrete categories.]
    ) <MechanicalEnergyGrouped>],
  [#figure(image("Graphics/BBShield Overlaid Mechanical Energy.png"),
  caption:[Bouncing ball shield overlaid with @MechanicalEnergyGrouped.]
  ) <BBShieldOverlaidMechanicalEnergy>],
  [#figure(image("Graphics/Possible to Hit.png"),
  caption:[Overlay to @MechanicalEnergyGrouped showing states where it is possible to hit the ball.]
  ) <PossibleToHit>]
)

@BBShieldOverlaidMechanicalEnergy shows the shield for the bouncing ball overlaid with this mechanical energy plot.
It is not very scientific, but I made sure the axis matched a visualisation of the safety strategy for the bouncing ball, and then I used a graphics editing program to overlay one with the other. Thus, the visualisation is not guaranteed to be pixel perfect, but it is very nearly correct.

@PossibleToHit shows the area where it is possible to hit the ball. Note that it is not known whether the ball will be hit simply from observing $E_mek$.

#pagebreak()

= Appendix

#figure(image("Graphics/BB Queries.png"), caption:[Learning queries.])<BBQueries>