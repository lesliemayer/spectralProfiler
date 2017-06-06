class DraggableAnnotation:
    """
    Adapted from: http://matplotlib.org/users/event_handling.html#draggable-rectangle-exercise
    """
    lock = None
    def __init__(self, annotation):
        self.ann = annotation
        self.press = None
        self.text = self.ann.get_text()

    def connect(self):
        self.cidpress = self.ann.figure.canvas.mpl_connect('button_press_event', self.onPress)
        self.cidrelease = self.ann.figure.canvas.mpl_connect('button_release_event', self.onRelease)
        self.cidmotion = self.ann.figure.canvas.mpl_connect('motion_notify_event', self.onMotion)

    def onPress(self, event):
        if event.inaxes != self.ann.axes:return

        contains, attrd = self.ann.contains(event)
        if not contains: return
        x0, y0 = self.ann.xy
        self.press = x0, y0, event.xdata, event.ydata
        DraggableAnnotation.lock = self

        #self.ann.set_text(r'\underline\{%f\}' % self.text)

        canvas = self.ann.figure.canvas
        axes = self.ann.axes
        self.ann.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.ann.axes.bbox)
        #print self.background
        axes.draw_artist(self.ann)

        canvas.blit(axes.bbox)

    def onRelease(self, event):
        #print "Releasing"
        #print self.ann.xy
        if DraggableAnnotation.lock is not self: return
        self.press = None
        DraggableAnnotation.lock = None

        self.ann.set_animated(False)
        self.background = None

        #Update the annotation position
        #self.ann.textcoords = (str(self.newx), str(self.newy))

        self.ann.figure.canvas.draw()

    def onMotion(self, event):
        if DraggableAnnotation.lock is not self: return
        if event.inaxes != self.ann.axes:return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.newx = x0 + dx
        self.newy = y0 + dy
        self.ann.xytext = (self.newx, self.newy)
        #print self.ann.xy
        canvas = self.ann.figure.canvas
        axes = self.ann.axes
        canvas.restore_region(self.background)
        axes.draw_artist(self.ann)

        canvas.blit(axes.bbox)

