graph [
  directed 1
  sustainability 1.0
  name "SCC1"
  node [
    id 0
    label "11"
    A 0
    C 1
    B 0
  ]
  node [
    id 1
    label "10"
    A 1
    C 0
    B 1
  ]
  node [
    id 2
    label "13"
    A 0
    C 1
    B 1
  ]
  node [
    id 3
    label "12"
    A 1
    C 1
    B 0
  ]
  node [
    id 4
    label "14"
    A 1
    C 1
    B 1
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
    internalweight 0.4
    weight 0.4
  ]
  edge [
    source 1
    target 2
    internalweight 0.6
    weight 0.6
  ]
  edge [
    source 2
    target 3
    internalweight 0.3
    weight 0.3
  ]
  edge [
    source 2
    target 4
    internalweight 0.7
    weight 0.7
  ]
  edge [
    source 3
    target 1
    internalweight 1.0
    weight 1.0
  ]
  edge [
    source 4
    target 2
    internalweight 1.0
    weight 1.0
  ]
]
