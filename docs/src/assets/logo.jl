using Luxor

Drawing(400, 400, "logo.svg")
origin()
background("white")
sethue("black")
squircle(Point(0,0), 150, 150, rt=0.3)
strokepath()
hoffset = 6.25
f(t) = Point(t-hoffset, cos(t)+3)
g(t) = Point(t-hoffset, cos(t+π/2))
h(t) = Point(t-hoffset, cos(t+π)-3)
setline(15)
fontsize(35)
@layer begin
    setopacity(0.4)
    sethue(Luxor.julia_red)
    poly(20f.(range(0, 4π, length=160)), :stroke)
    sethue(Luxor.julia_green)
    poly(20g.(range(0, 4π, length=160)), :stroke)
    sethue(Luxor.julia_purple)
    poly(20h.(range(0, 4π, length=160)), :stroke)
end
fontsize(135)
fontface("Helvetica")
text("WW", halign=:center, valign=:middle)
finish()
preview()
