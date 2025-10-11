from graphviz import Digraph

# More exhaustive version including return to S from intermediates
dot = Digraph(format='svg')
dot.attr(rankdir='LR', size='14,7')
dot.attr('node', shape='box', style='filled', fillcolor='lightgrey', fontsize='12')

# Nodes
dot.node('S', 'S\n(Single-stranded DNA with RAD51)', fillcolor='lightgrey')
dot.node('HM', 'Homologous Complexes (HM)', fillcolor='lightblue')
dot.node('HT', 'Heterologous Complexes (HT)', fillcolor='lightcoral')
dot.node('HMext', 'Extension Path\n(HM8 → HM9 → ... → HM192)', fillcolor='lightblue')
dot.node('HText', 'Extension Path\n(HT → HT9 → ... → HT192)', fillcolor='lightcoral')
dot.node('DHM', 'D-loop (Homologous)', fillcolor='lightblue')
dot.node('DHT', 'D-loop (Heterologous)', fillcolor='lightcoral')
dot.node('R', 'Final Recombinant State (R)', fillcolor='gold')
dot.node('ReturnS', 'Return to S', fillcolor='white')

# Edges from S to initial complexes
dot.edge('S', 'HM', label='kon × f\n(association)', fontsize='10')
dot.edge('S', 'HT', label='kon × (1 - f)\n(association)', fontsize='10')

# Extensions
dot.edge('HM', 'HMext', label='Stepwise extension\n(kext)', fontsize='10')
dot.edge('HT', 'HText', label='Stepwise extension\n(kext × ktol)', fontsize='10')

# D-loop formation
dot.edge('HMext', 'DHM', label='D-loop formation\n(kdloop, position-dependent)', fontsize='10')
dot.edge('HText', 'DHT', label='D-loop formation\n(kdloop × ktol)', fontsize='10')

# Recombination
dot.edge('DHM', 'R', label='Recombination\n(kre)', fontsize='10')
dot.edge('DHT', 'R', label='Recombination\n(kre)', fontsize='10')

# Dissociation from D-loops
dot.edge('DHM', 'ReturnS', label='koff2', style='dashed', fontsize='10')
dot.edge('DHT', 'ReturnS', label='koff2', style='dashed', fontsize='10')

# Dissociation from intermediates
dot.edge('HMext', 'ReturnS', label='Intermediate dissociation\n(koff1 / 1.4^n)', style='dashed', fontsize='10')
dot.edge('HText', 'ReturnS', label='Intermediate dissociation\n(koff1 / 1.4^n)', style='dashed', fontsize='10')

# Output path
svg_full_path = "gillespie_model"
dot.render(svg_full_path, cleanup=False)
svg_full_path

