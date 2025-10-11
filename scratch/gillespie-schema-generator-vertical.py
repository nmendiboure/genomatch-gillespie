from graphviz import Digraph

# Airier layout with more spacing between nodes, still in 4:3 ratio
dot = Digraph(format='svg')
dot.attr(rankdir='TB', size='8,6')  # Top-to-Bottom layout, 4:3 ratio
dot.attr('graph', nodesep='1.0', ranksep='1.2')  # Increase spacing
dot.attr('node', shape='box', style='filled', fillcolor='lightgrey', fontsize='14', margin='0.3,0.2')

# Nodes
dot.node('S', 'S\n(Single-stranded DNA with RAD51)', fillcolor='lightgrey')
dot.node('HM', 'Homologous Complexes (HM)', fillcolor='lightblue')
dot.node('HT', 'Heterologous Complexes (HT)', fillcolor='lightcoral')
dot.node('HMext', 'Extension Path\n(HM8 → HM9 → ... → HM195)', fillcolor='lightblue')
dot.node('HText', 'Extension Path\n(HT8 → HT9 → ... → HT195)', fillcolor='lightcoral')
dot.node('DHM', 'D-loop\n(Homologous)', fillcolor='lightblue')
dot.node('DHT', 'D-loop\n(Heterologous)', fillcolor='lightcoral')
dot.node('R', 'Final Recombinant\nState (R)', fillcolor='gold')
dot.node('ReturnS', 'Return to S', fillcolor='white')

# Edges with spacing-friendly formatting
dot.edge('S', 'HM', label='kon × f\n(association)', fontsize='14')
dot.edge('S', 'HT', label='kon × (1 - f)\n(association)', fontsize='14')

dot.edge('HM', 'HMext', label='Extension\n(kext)', fontsize='14')
dot.edge('HT', 'HText', label='Extension\n(kext × ktol)', fontsize='14')

dot.edge('HMext', 'DHM', label='D-loop formation\n(kdloop)', fontsize='14')
dot.edge('HText', 'DHT', label='D-loop formation\n(kdloop × ktol)', fontsize='14')

dot.edge('DHM', 'R', label='Recombination\n(kre)', fontsize='14')
dot.edge('DHT', 'R', label='Recombination\n(kre)', fontsize='14')

dot.edge('DHM', 'ReturnS', label='koff2', style='dashed', fontsize='14')
dot.edge('DHT', 'ReturnS', label='koff2', style='dashed', fontsize='14')

dot.edge('HMext', 'ReturnS', label='Intermediate dissociation\n(koff1 / 1.4^n)', style='dashed', fontsize='14')
dot.edge('HText', 'ReturnS', label='Intermediate dissociation\n(koff1 / 1.4^n)', style='dashed', fontsize='14')

# Output SVG path
svg_aerated_path = "gillespie_model_4x3"
dot.render(svg_aerated_path, cleanup=False)
svg_aerated_path

