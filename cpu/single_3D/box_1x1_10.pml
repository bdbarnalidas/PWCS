set float_labels,on
set label_color,blue
pseudoatom a1, pos=[-0.500000, 10.500000, -0.500000]
pseudoatom a2, pos=[10.500000, 10.500000, -0.500000]
distance d1, /a1/////1, /a2/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d1
hide labels, d1
pseudoatom a3, pos=[10.500000, 10.500000, 10.500000]
distance d2, /a2/////1, /a3/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d2
hide labels, d2
pseudoatom a4, pos=[-0.500000, 10.500000, 10.500000]
distance d3, /a3/////1, /a4/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d3
hide labels, d3
distance d4, /a1/////1, /a4/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d4
hide labels, d4
pseudoatom a5, pos=[-0.500000, -0.500000, -0.500000]
pseudoatom a6, pos=[10.500000, -0.500000, -0.500000]
distance d5, /a5/////1, /a6/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d5
hide labels, d5
pseudoatom a7, pos=[10.500000, -0.500000, 10.500000]
distance d6, /a6/////1, /a7/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d6
hide labels, d6
pseudoatom a8, pos=[-0.500000, -0.500000, 10.500000]
distance d7, /a7/////1, /a8/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d7
hide labels, d7
distance d8, /a5/////1, /a8/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d8
hide labels, d8
distance d9, /a1/////1, /a5/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d9
hide labels, d9
distance d10, /a2/////1, /a6/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d10
hide labels, d10
distance d11, /a4/////1, /a8/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d11
hide labels, d11
distance d12, /a3/////1, /a7/////1
set dash_gap, 0
set dash_color, 0xadd8e6, d12
hide labels, d12
show nb_spheres
color 0xadd8e6
