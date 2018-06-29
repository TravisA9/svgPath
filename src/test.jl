using Gtk, Graphics, Cairo, Colors
using svgPath

include("samples.jl")
path = parsePathData(Computer)
#path = ScalePath(path, 2)



c = @GtkCanvas()
win = GtkWindow(c, "Canvas")
@guarded draw(c) do widget
    ctx = getgc(c)
    h = height(c)
    w = width(c)
    rectangle(ctx, 50,50, 50, 50)
    set_source_rgb(ctx, 1, 1, 1)
    fill(ctx)

    set_line_width(ctx, 1)
    set_source_rgb(ctx, 0, 0, 0)
    # stroke(ctx)
    drawPath(ctx, path)
    #println(get_current_point(ctx))
    #println(copy_path_flat(ctx))

    stroke(ctx)
end
show(c)
