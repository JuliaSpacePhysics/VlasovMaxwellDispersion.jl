#import "@preview/cetz:0.4.2"

#cetz.canvas(length: 1cm, {
  import cetz.draw: *
  let d = (2.1, -0.75)
  let ray(a, r) = (2.1 + r * calc.cos(a), -0.75 + r * calc.sin(a))
  line(d, ray(-18deg, 2.9), ray(18deg, 3.05), ray(55deg, 2.9), close: true, fill: rgb("fdf3df"), stroke: none)
  line(d, ray(-18deg, 2.9), stroke: (dash: "dashed", thickness: 0.5pt))
  line(d, ray(55deg, 2.9), stroke: (dash: "dashed", thickness: 0.5pt))
  content(ray(40deg, 2.45), $S(omega)$)
  line((-0.5, 0), (5.3, 0), stroke: 0.6pt)
  content((5.5, 0.32), $"Re" y$)
  line((0, -1.6), (0, 1.6), stroke: 0.6pt)
  content((0.55, 1.62), $"Im" y$)
  line((1.6, 0), (4.9, 0), mark: (end: "stealth"), stroke: (paint: gray, thickness: 1.3pt))
  circle((1.6, 0), radius: 0.07, stroke: 0.6pt, fill: white)
  content((3.5, 0.2), text(fill: gray.darken(20%), size: 8pt)[germ ray])
  bezier((1.6, 0), d, (1.5, -0.6), mark: (end: "stealth", fill: black), stroke: 0.7pt)
  circle(d, radius: 0.07, fill: rgb("c0392b"), stroke: none)
  content((1.5, -1.15), $y_0(omega)$)
  line(d, ray(20deg, 2.8), mark: (end: "stealth", fill: rgb("2457a6")), stroke: (paint: rgb("2457a6")))
  content((3.85, -0.5), text(fill: rgb("2457a6"))[$Gamma(omega)$])
})
