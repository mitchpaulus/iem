#import "@preview/lilaq:0.5.0" as lq
#import "@preview/cetz:0.4.2"

// ---------------------------------------------------------------------------
// Style
// ---------------------------------------------------------------------------
#let ink = rgb("#0b0b0b")
#let ink2 = rgb("#52514e")
#let muted = rgb("#8a8984")
#let c1 = rgb("#2a78d6") // blue   - DFW / station pressure
#let c2 = rgb("#eb6834") // orange - DEN / altimeter setting
#let c3 = rgb("#1baf7a") // aqua

#set document(title: "Altimeter Setting vs. Station Pressure")
#set page(paper: "us-letter", margin: (x: 1in, y: 0.9in), numbering: "1")
#set text(size: 10.5pt, fill: ink)
#set par(justify: true)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show figure.caption: set text(size: 9pt, fill: ink2)
#show: lq.set-diagram(width: 12cm, height: 6cm)
#show: lq.set-legend(stroke: 0.5pt + muted.lighten(40%))
#show: lq.set-grid(stroke: 0.4pt + muted.lighten(65%))

// ---------------------------------------------------------------------------
// Physics
// ---------------------------------------------------------------------------
#let T0 = 518.67          // °R, standard sea-level temperature (59 °F)
#let lapse = 0.0035662    // °R/ft, standard lapse rate (3.566 °F per 1000 ft)
#let Rd = 53.353          // ft·lbf/(lbm·°R), dry air
#let expo-derived = 1 / (Rd * lapse)   // (g/gc)/(R alpha), g/gc = 1 lbf/lbm
#let coef-derived = lapse / T0         // 1/ft
// ASHRAE Fundamentals Ch. 1 Eq. 3 (IP) published constants, used everywhere else.
#let expo = 5.2559
#let coef = 6.8754e-6
#let psi-per-inhg = 0.491154   // psi per inHg

// Station pressure from altimeter setting, standard-atmosphere (ASHRAE) form.
// Pressure in any unit, h in ft.
#let p-ashrae(alti, h) = alti * calc.pow(1 - coef * h, expo)
// Station pressure from altimeter setting, NWS form converted to inHg and ft.
#let nws-n = 0.190284
#let nws-k-si = calc.pow(1013.25, nws-n) * 0.0065 / 288        // hPa^n per m
#let nws-k = nws-k-si * 0.3048 / calc.pow(33.8639, nws-n)      // inHg^n per ft
#let nws-d = 0.3 / 33.8639                                     // inHg
#let p-nws(alti, h) = calc.pow(calc.pow(alti, nws-n) - nws-k * h, 1 / nws-n) + nws-d
// Saturation pressure over liquid water, psia (ASHRAE Hyland-Wexler, IP). t in °F.
#let pws(t) = {
  let T = t + 459.67
  calc.exp(
    -1.0440397e4 / T - 1.1294650e1 - 2.7022355e-2 * T + 1.2890360e-5 * T * T
      - 2.4780681e-9 * T * T * T + 6.5459673 * calc.ln(T),
  )
}
#let W(pv, p) = 0.621945 * pv / (p - pv)

#let fmt(x, d) = str(calc.round(x, digits: d))

// Relative errors from using alti as absolute pressure (alti = 29.921 inHg, 60 °F dew point)
#let a-std = 29.921
#let perr-at(h) = (a-std / p-ashrae(a-std, h) - 1) * 100
#let werr-at(h) = {
  let pv = pws(60.0)
  let pa = a-std * psi-per-inhg
  (1 - W(pv, pa) / W(pv, p-ashrae(pa, h))) * 100
}

// ---------------------------------------------------------------------------
// Data: IEM ASOS, September 2026, UTC, routine (report_type=3) observations
// ---------------------------------------------------------------------------
#let H = (DFW: 182.0 * 3.28084, DEN: 1656.0 * 3.28084) // ft (IEM metadata: 182 m, 1656 m)
#let raw = csv("data/dfw_den_2026-09.csv").slice(1)
#let obs(stn) = {
  raw.filter(r => r.at(0) == stn and r.at(4) != "M").map(r => {
    let v = r.at(1)
    (
      day: int(v.slice(8, 10)) + int(v.slice(11, 13)) / 24 + int(v.slice(14, 16)) / 1440,
      tmpf: if r.at(2) == "M" { none } else { float(r.at(2)) },
      dwpf: if r.at(3) == "M" { none } else { float(r.at(3)) },
      alti: float(r.at(4)),                                             // inHg
      mslp: if r.at(5) == "M" { none } else { float(r.at(5)) / 33.8639 }, // mb -> inHg
    )
  })
}
#let dfw = obs("DFW")
#let den = obs("DEN")

// Worked example: first DFW observation
#let ex = dfw.at(0)
#let ex-pstn = p-ashrae(ex.alti, H.DFW)          // inHg
#let ex-pv = pws(ex.dwpf)                         // psia
#let ex-W-true = W(ex-pv, ex-pstn * psi-per-inhg)
#let ex-W-alti = W(ex-pv, ex.alti * psi-per-inhg)

// ---------------------------------------------------------------------------
#align(center)[
  #text(size: 17pt, weight: "bold")[Altimeter Setting vs. Station Pressure]
  #v(-0.4em)
  #text(size: 11pt, fill: ink2)[
    Getting absolute air pressure for psychrometrics from IEM ASOS `alti`
  ]
]
#v(0.6em)

#block(
  fill: c1.lighten(90%), inset: 10pt, radius: 3pt, width: 100%,
)[
  *Summary.* The ASOS `alti` field is _not_ the air pressure at the station. It is
  the station pressure adjusted to sea level with a fixed, idealized atmosphere.
  Because that atmosphere is fixed, the adjustment can be undone exactly:
  $ p_"stn" = p_"alt" (1 - #fmt(coef * 1e6, 4) times 10^(-6) H)^(#fmt(expo, 4)) $
  with $H$ the station elevation in feet. Feeding `alti` straight into
  psychrometric equations overstates pressure by about #fmt(coef * expo * 1e5, 1) % per 1000 ft
  near sea level. At Denver that makes the humidity ratio about #fmt(werr-at(H.DEN), 0) % too low.
]

= What the altimeter setting is

An ASOS station measures the true (absolute) air pressure at the elevation of its
barometer, the *station pressure* $p_"stn"$. Pilots need something else: a number
to dial into an aircraft altimeter so that, on the runway, it reads the field
elevation. That number is the *altimeter setting* $p_"alt"$, reported in METARs
as `A2991` (29.91 inHg) and in IEM data as `alti`.

The altimeter is just a barometer with a height scale. Its scale assumes the
*ICAO standard atmosphere*. The altimeter setting is defined as the sea-level
pressure of the standard atmosphere that passes through the measured
$p_"stn"$ at the station's elevation $H$ (@fig-schematic).

#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *
    // sea level
    line((0, 0), (13, 0), stroke: 1pt + c1.darken(20%))
    content((0, -0.15), anchor: "north-west", text(8pt, fill: c1.darken(20%))[sea level, $z = 0$])
    // terrain
    line(
      (0, 0), (1.5, 0.3), (3, 1.0), (4.5, 2.3), (6, 3), (9, 3), (10.5, 2.4), (13, 2.0),
      fill: muted.lighten(80%), stroke: 0.8pt + muted,
    )
    // station
    rect((7.1, 3), (7.9, 3.9), fill: white, stroke: 0.8pt + ink)
    content((7.5, 3.45), text(7pt)[ASOS])
    line((7.5, 3.9), (7.5, 4.6), stroke: 0.8pt + ink)
    content((7.5, 4.65), anchor: "south", text(8.5pt)[measures $p_"stn"$ at $z = H$])
    // elevation dimension
    line((11.6, 0), (11.6, 3), mark: (start: "stealth", end: "stealth", fill: ink2), stroke: 0.6pt + ink2)
    content((11.75, 1.5), anchor: "west", text(8.5pt, fill: ink2)[$H$])
    line((7.9, 3), (11.8, 3), stroke: (paint: ink2, thickness: 0.4pt, dash: "dotted"))
    // imaginary standard-atmosphere column
    rect((6.6, 0), (8.4, 3), fill: c2.lighten(85%), stroke: (paint: c2, thickness: 0.8pt, dash: "dashed"))
    line((7.5, 2.85), (7.5, 0.25), mark: (end: "stealth", fill: c2), stroke: 1.2pt + c2)
    content((6.5, 1.5), anchor: "east", align(right, text(8.5pt)[imaginary column of\ *standard atmosphere*\ $T = 59 degree"F" - 3.566 degree"F"\/1000 "ft" dot z$]))
    content((7.5, -0.15), anchor: "north", text(8.5pt)[pressure at the bottom $= p_"alt"$ (`alti`)])
    // real atmosphere note
    content((2.0, 3.7), anchor: "west", align(left, text(8pt, fill: ink2)[the real air temperature\ is *not* used]))
  }),
  caption: [
    The altimeter setting extends the station down to sea level through an
    imaginary column of standard atmosphere. The real air temperature never
    enters, which is what makes the step exactly reversible.
  ],
) <fig-schematic>

= The math

== Pressure in the standard atmosphere

The standard atmosphere has sea-level temperature $T_0 = 518.67 degree"R"$ (59 °F)
and a constant lapse rate $alpha = 0.0035662 degree"R"\/"ft"$ in the troposphere:
$ T(z) = T_0 - alpha z. $
Combining hydrostatic balance with the ideal-gas law for dry air,
$ (dif p) / (dif z) = - rho g / g_c = - p / (R_d T(z)) dot g / g_c, $
and integrating from sea level ($p = p_0$) up to height $z$ gives
$ p(z) = p_0 (1 - (alpha z) / T_0)^((g \/ g_c) \/ (R_d alpha)). $ <eq-std>
With standard gravity $g \/ g_c = 1$ lbf/lbm and $R_d = 53.353$ ft·lbf/(lbm·°R)
the constants are
$ alpha / T_0 = #fmt(coef-derived * 1e6, 4) times 10^(-6) "ft"^(-1), quad
  (g \/ g_c) / (R_d alpha) = #fmt(expo-derived, 4). $
@eq-std is the same equation ASHRAE _Fundamentals_ (Ch. 1, Eq. 3) uses for
standard pressure at altitude, with $p_0 = 14.696$ psia. ASHRAE publishes the
constants as $6.8754 times 10^(-6)$ and $5.2559$ (the same to rounding), and the
rest of this document uses those values.

== Defining and inverting the altimeter setting

The altimeter setting is the $p_0$ for which @eq-std passes through the
measured station pressure at the station elevation:
$ p_"stn" = p_"alt" (1 - (alpha H) / T_0)^((g \/ g_c) \/ (R_d alpha))
  quad <==> quad
  p_"alt" = p_"stn" (1 - (alpha H) / T_0)^(-(g \/ g_c) \/ (R_d alpha)). $ <eq-inv>
The factor depends only on $H$, so it is a fixed constant for each station.
Station pressure is `alti` times that constant (@fig-profile).

#figure(
  {
    let a-dfw = dfw.at(0).alti
    let a-den = den.at(0).alti
    let zs = lq.linspace(0, 6500, num: 60)
    let s-dfw = p-ashrae(a-dfw, H.DFW)
    let s-den = p-ashrae(a-den, H.DEN)
    lq.diagram(
      xlabel: [Pressure (inHg)],
      ylabel: [Elevation $z$ (ft)],
      xlim: (23, 31),
      ylim: (-250, 6500),
      yaxis: (exponent: 0),
      legend: (position: top + right),
      lq.plot(zs.map(z => p-ashrae(a-dfw, z)), zs, mark: none, color: c1,
        stroke: 1.5pt + c1, label: [DFW standard-atmosphere curve]),
      lq.plot(zs.map(z => p-ashrae(a-den, z)), zs, mark: none, color: c2,
        stroke: (paint: c2, thickness: 1.5pt, dash: "dashed"), label: [DEN standard-atmosphere curve]),
      lq.hlines(0, stroke: 0.6pt + ink2),
      lq.hlines(H.DFW, H.DEN, stroke: (paint: muted, thickness: 0.5pt, dash: "dotted")),
      lq.scatter((s-dfw, s-den), (H.DFW, H.DEN), color: (c1, c2), size: 7pt, mark: "o"),
      lq.scatter((a-dfw, a-den), (0, 0), color: (c1, c2), size: 7pt, mark: "s"),
      lq.place(s-dfw - 0.24, H.DFW, align: right + horizon,
        text(8pt)[DFW $p_"stn"$ = #fmt(s-dfw, 2)]),
      lq.place(s-den - 0.24, H.DEN, align: right + horizon,
        text(8pt)[DEN $p_"stn"$ = #fmt(s-den, 2)]),
      lq.place(28.35, 100, align: right + bottom,
        text(8pt)[`alti`: DFW #fmt(a-dfw, 2), DEN #fmt(a-den, 2)]),
    )
  },
  caption: [
    @eq-std drawn through the first September 2026 observation at DFW (#fmt(H.DFW, 0) ft) and
    Denver (DEN, #fmt(H.DEN, 0) ft). Circles: station pressure at the station elevation.
    Squares: where each curve crosses sea level, which is the reported `alti`. Both
    stations report about the same `alti`, but Denver's station pressure is about
    #fmt((1 - p-ashrae(den.at(0).alti, H.DEN) / p-ashrae(dfw.at(0).alti, H.DFW)) * 100, 0) % lower.
  ]
) <fig-profile>

== The NWS formula

The NWS altimeter-setting formula is the same relation with
$n = 0.190284 approx 1 \/ 5.2553$. With pressure in inHg and $H$ in feet it reads
$ p_"alt"^n = (p_"stn" - d)^n + K H, $
so
$ p_"stn" = (p_"alt"^n - K H)^(1 \/ n) + d, quad
  K = #fmt(nws-k * 1e5, 4) times 10^(-5), quad d = #fmt(nws-d, 4) "inHg". $ <eq-nws>
Apart from slightly rounded constants and the small fixed offset $d$, this is @eq-inv. @fig-formulas shows the two forms differ by less than
0.009 inHg below 10,000 ft, mostly from the offset. That is under one reporting
step of `alti` (0.01 inHg) and about 0.03 % of the pressure, so either form is fine.

#figure(
  {
    let hs = lq.linspace(0, 10000, num: 61)
    let a = a-std
    lq.diagram(
      height: 4.5cm,
      xlabel: [Station elevation $H$ (ft)],
      ylabel: [NWS − ASHRAE (inHg)],
      ylim: (0, 0.015),
      xaxis: (exponent: 0),
      legend: (position: top + right),
      lq.fill-between(hs, hs.map(_ => 0), y2: hs.map(_ => 0.01),
        fill: muted.lighten(80%), stroke: none, label: [one `alti` step (0.01 inHg)]),
      lq.plot(hs, hs.map(h => p-nws(a, h) - p-ashrae(a, h)), mark: none,
        stroke: 1.5pt + c1, label: [difference in computed $p_"stn"$, $p_"alt"$ = 29.921 inHg]),
    )
  },
  caption: [
    Difference between @eq-nws and @eq-inv. It is smaller than one reporting
    step of `alti` everywhere on the plot.
  ],
) <fig-formulas>

= Why it matters for psychrometrics

The humidity ratio is
$ W = 0.621945 p_v / (p - p_v), $
in lb#sub[w]/lb#sub[da], where $p_v$ is the water vapor partial pressure (the
saturation pressure at the dew point) and $p$ is the *absolute* total pressure,
both in psia. Using $p_"alt"$ in place of $p_"stn"$ gives a
relative error of
$ W_"alt" / W_"true" - 1 = (p_"stn" - p_v) / (p_"alt" - p_v) - 1 approx -(p_"alt" - p_"stn") / p_"alt". $
So the humidity ratio is understated by about the same percentage that
pressure is overstated. Enthalpy, specific volume, and air density inherit
the error. Wet-bulb temperature is affected less directly.

*Worked example: DFW, 1 Sep 2026 00:53 UTC.* Dry bulb
#fmt(ex.tmpf, 0) °F, dew point #fmt(ex.dwpf, 0) °F, `alti` = #fmt(ex.alti, 2) inHg, $H$ = #fmt(H.DFW, 0) ft.

#table(
  columns: (1fr, auto, auto),
  align: (left, right, right),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  table.hline(stroke: 0.6pt),
  table.header([], [using `alti`], [using $p_"stn"$]),
  table.hline(stroke: 0.4pt),
  [Pressure (inHg)], [#fmt(ex.alti, 2)], [#fmt(ex-pstn, 2)],
  [Pressure (psia)], [#fmt(ex.alti * psi-per-inhg, 3)], [#fmt(ex-pstn * psi-per-inhg, 3)],
  [$p_v$ at dew point (psia)], [#fmt(ex-pv, 4)], [#fmt(ex-pv, 4)],
  [Humidity ratio (lb/lb)], [#fmt(ex-W-alti, 5)], [#fmt(ex-W-true, 5)],
  [Humidity ratio (gr/lb)], [#fmt(ex-W-alti * 7000, 1)], [#fmt(ex-W-true * 7000, 1)],
  table.hline(stroke: 0.6pt),
)
Using `alti` gives a humidity ratio #fmt((1 - ex-W-alti / ex-W-true) * 100, 1) % too low at DFW.
@fig-error shows how the error grows with elevation.

#figure(
  {
    let hs = lq.linspace(0, 8000, num: 51)
    let perr = perr-at
    let werr = werr-at
    lq.diagram(
      xlabel: [Station elevation $H$ (ft)],
      ylabel: [Error from using `alti` (%)],
      xlim: (0, 8000),
      ylim: (-3, 40),
      xaxis: (exponent: 0),
      legend: (position: top + left),
      lq.plot(hs, hs.map(perr), mark: none, stroke: 1.5pt + c1,
        label: [pressure overstated]),
      lq.plot(hs, hs.map(werr), mark: none,
        stroke: (paint: c2, thickness: 1.5pt, dash: "dashed"),
        label: [humidity ratio understated (60 °F dew point)]),
      lq.vlines(H.DFW, H.DEN, stroke: (paint: muted, thickness: 0.5pt, dash: "dotted")),
      lq.scatter((H.DFW, H.DEN), (perr(H.DFW), perr(H.DEN)), color: c1, size: 6pt),
      lq.scatter((H.DFW, H.DEN), (werr(H.DFW), werr(H.DEN)), color: c2, size: 6pt, mark: "s"),
      lq.place(H.DFW + 400, -1.4, align: left + horizon,
        text(8pt)[DFW: $p$ +#fmt(perr(H.DFW), 1) %, $W$ −#fmt(werr(H.DFW), 1) %]),
      lq.place(H.DEN + 200, werr(H.DEN) - 1, align: left + top,
        text(8pt)[DEN: $p$ +#fmt(perr(H.DEN), 1) %, $W$ −#fmt(werr(H.DEN), 1) %]),
    )
  },
  caption: [
    Relative error from treating `alti` as absolute pressure. The pressure
    error is $(1 - #fmt(coef * 1e6, 4) times 10^(-6) H)^(-#fmt(expo, 4)) - 1$ with $H$ in ft. The humidity-ratio
    error is about the same size, slightly smaller, and nearly independent of dew point.
  ],
) <fig-error>

= What it looks like in real data

@fig-series compares `alti` and the station pressure from @eq-inv over
September 2026. Both move with the weather, but the gap is the same fixed
ratio throughout. At Denver the gap is several times larger than the
month's whole weather range.

#figure(
  {
    let panel(rows, h, name, ylim) = lq.diagram(
      width: 12cm, height: 3.6cm,
      title: [#name ($H$ = #fmt(h, 0) ft)],
      ylabel: [psia],
      xlim: (1, 31),
      ylim: ylim,
      legend: none,
      lq.plot(rows.map(r => r.day), rows.map(r => r.alti * psi-per-inhg), mark: none,
        stroke: (paint: c2, thickness: 1pt, dash: "dashed")),
      lq.plot(rows.map(r => r.day), rows.map(r => p-ashrae(r.alti, h) * psi-per-inhg), mark: none,
        stroke: 1pt + c1),
      lq.place(1.3, rows.at(0).alti * psi-per-inhg, align: left + bottom,
        box(inset: (bottom: 5pt), text(8pt)[`alti` (as if absolute)])),
      lq.place(1.3, p-ashrae(rows.at(0).alti, h) * psi-per-inhg, align: left + top,
        box(inset: (top: 5pt), text(8pt)[$p_"stn"$ (true absolute pressure)])),
    )
    stack(spacing: 0.6em,
      panel(dfw, H.DFW, [DFW], (14.1, 14.9)),
      panel(den, H.DEN, [DEN], (11.5, 15.5)),
      align(center, text(9pt)[Day of September 2026 (UTC)]),
    )
  },
  caption: [
    Hourly `alti` (dashed) and station pressure (solid) in psia.
  ],
) <fig-series>

= Why not use `mslp`?

IEM also offers `mslp`, the sea-level pressure from the METAR remarks (`SLP112`
→ 1011.2 mb = 29.86 inHg; IEM reports it in mb). It sounds like the more
"scientific" choice, but it is worse for this purpose:

- It is reduced to sea level with the *actual* air temperature (the mean of the
  current and 12-hour-old temperatures), plus empirical corrections at high
  elevation. Undoing it requires those same inputs and corrections.
- It is optional in METARs, so it is missing more often than `alti`.

@fig-mslp shows the temperature dependence directly. At DFW, `mslp` and
`alti` differ by a nearly constant 0.03 inHg. At Denver, the difference swings by
0.1–0.3 inHg as the temperature changes, even though both describe the same
station pressure.

#figure(
  {
    let pts(rows) = rows.filter(r => r.mslp != none and r.tmpf != none)
    let d1 = pts(dfw)
    let d2 = pts(den)
    lq.diagram(
      xlabel: [Air temperature (°F)],
      ylabel: [`mslp` − `alti` (inHg)],
      legend: (position: top + left),
      lq.hlines(0, stroke: 0.6pt + ink2),
      lq.scatter(d1.map(r => r.tmpf), d1.map(r => r.mslp - r.alti),
        color: c1, size: 3.5pt, mark: "o", label: [DFW (#fmt(H.DFW, 0) ft)]),
      lq.scatter(d2.map(r => r.tmpf), d2.map(r => r.mslp - r.alti),
        color: c2, size: 3.5pt, mark: "s", label: [DEN (#fmt(H.DEN, 0) ft)]),
    )
  },
  caption: [
    `mslp` − `alti` vs. temperature, hourly, September 2026. The altimeter setting
    has no temperature term, but sea-level pressure does. At high elevation the
    two disagree by an amount that depends on the weather.
  ],
) <fig-mslp>

= Recipe

+ Get `alti` in inHg, as reported.
+ Get the station elevation $H$ in feet. IEM publishes it in meters in the network
  metadata, e.g. `geojson/network/TX_ASOS.geojson`. Multiply by 3.28084
  (DFW: #fmt(H.DFW, 0) ft, DEN: #fmt(H.DEN, 0) ft).
+ $p_"stn" ["inHg"] = p_"alt" (1 - #fmt(coef * 1e6, 4) times 10^(-6) H)^(#fmt(expo, 4))$.
+ For psychrometrics, convert: psia = 0.491154 × inHg.

#v(1em)
#text(size: 8.5pt, fill: ink2)[
  Data: Iowa Environmental Mesonet ASOS archive (`/cgi-bin/request/asos.py`),
  stations DFW and DEN, 2026-09-01 to 2026-10-01 UTC, routine observations.
  Vapor pressure is the saturation pressure over liquid water at the dew point,
  from the Hyland–Wexler equation in IP units (ASHRAE _Fundamentals_, Ch. 1, Eq. 6),
  $T$ in °R, $p_"ws"$ in psia.
]
