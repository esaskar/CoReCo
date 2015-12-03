def plot_format(plotfn, format):
    if format == "png":
        return "png(file=\"%s.png\",width=1024,height=1024,pointsize=24)\n" % (plotfn)
    elif format == "pdf":
        return "pdf(file=\"%s.pdf\")\n" % (plotfn)
    else:
        raise Exception("Unknown format: %s" % (format))

