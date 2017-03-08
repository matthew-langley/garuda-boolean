graph [
  directed 1
  name "SimpleGraphWithSCCs"
  node [
    id 0
    label 0
	A 0
	B 0
	C 0
  ]
  node [
    id 1
    label 1
	A 0
	B 0
	C 0
  ]
  node [
    id 2
    label 2
	A 1
	B 0
	C 0
  ]
  node [
    id 3
    label 3
	A 0
	B 1
	C 0
  ]
  node [
    id 4
    label 10
	A 1
	B 1
	C 0
  ]
  node [
    id 5
    label 11
	A 0
	B 0
	C 1
  ]
  node [
    id 6
    label 12
	A 1
	B 0
	C 1
  ]
  node [
    id 7
    label 13
	A 0
	B 1
	C 1
  ]
  node [
    id 8
    label 14
	A 1
	B 1
	C 1
  ]
  edge [
    source 0
    target 1
    weight 0.8
  ]
  edge [
    source 0
    target 4
    weight 0.2
  ]
  edge [
    source 1
    target 2
    weight 1.0
  ]
  edge [
    source 2
    target 3
    weight 1.0
  ]
  edge [
    source 3
    target 0
    weight 1.0
  ]
  edge [
    source 4
    target 5
    weight 0.4
  ]
  edge [
    source 4
    target 7
    weight 0.6
  ]
  edge [
    source 5
    target 6
    weight 1.0
  ]
  edge [
    source 6
    target 4
    weight 1.0
  ]
  edge [
    source 7
    target 6
    weight 0.3
  ]
  edge [
    source 7
    target 8
    weight 0.7
  ]
  edge [
    source 8
    target 7
    weight 1.0
  ]
]
