from ete3 import Tree


tree = Tree("protein_tree.nwk")

print("Original tree:")
print(tree)


camel_virus = None
for node in tree.iter_leaves():
    if "MN514967" in node.name:
        camel_virus = node
        break

if camel_virus:
    print(f"\nFound outgroup: {camel_virus.name}")
    
    tree.set_outgroup(camel_virus)
    
    print("\nRooted tree:")
    print(tree)

    tree.write(outfile="rooted_protein_tree.nwk")
    print("\nRooted tree saved to: rooted_protein_tree.nwk")
else:
    print("Camel virus (MN514967) not found in the tree!")
    print("Available leaves:")
    for leaf in tree.iter_leaves():
        print(f"  - {leaf.name}")
