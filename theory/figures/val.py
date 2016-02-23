Uxmax = 0.0002099
Uxmin = -0.008879

Uymax = 0.03136
Uymin = -0.02856
DragMax = 474.6
DragMin = 408.3

LiftMax = 175.7
LiftMin = -180.3


Ux_mean = (Uxmax+Uxmin)*0.5
Ux_amp = (Uxmax-Uxmin)*0.5

Uy_mean = (Uymax+Uymin)*0.5
Uy_amp = (Uymax - Uymin)*0.5

Drag_mean = (DragMax + DragMin)*0.5
Drag_amp = (DragMax - DragMin)*0.5

Lift_mean = (LiftMax + LiftMin)*0.5
Lift_amp = (LiftMax - LiftMin)*0.5

print 'Ux = %.2f \pm %.2f' %(Ux_mean*1e3,Ux_amp*1e3)
print 'Uy = %.2f \pm %.2f'%(Uy_mean*1e3,Uy_amp*1e3)
print 'Drag = %.2f \pm %.2f'%(Drag_mean,Drag_amp)
print 'Lift = %.2f \pm %.2f'%(Lift_mean,Lift_amp)
