graph [
  directed 1
  sustainability 0.875
  name "SCC2"
  node [
    id 0
    label "1"
    A 0
    C 0
    B 0
  ]
  node [
    id 1
    label "0"
    A 0
    C 0
    B 0
  ]
  node [
    id 2
    label "3"
    A 0
    C 0
    B 1
  ]
  node [
    id 3
    label "2"
    A 1
    C 0
    B 0
  ]
  edge [
    source 0
    target 3
    internalweight 1.0
    weight 1.0
  ]
  edge [
    source 1
    target 0
    internalweight 1.0
    weight 0.8
  ]
  edge [
    source 2
    target 1
    internalweight 1.0
    weight 1.0
  ]
  edge [
    source 3
    target 2
    internalweight 1.0
    weight 1.0
  ]
]
