from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Define the reaction using SMARTS
reaction_smarts = "[CH3CH2O-].[Na+].[CH3CH2CH2CH2Br]>>[CH3CH2OCH2CH2CH2CH3].[Na+].[Br-]"
reaction = AllChem.ReactionFromSmarts(reaction_smarts)

# Define reactants as RDKit molecules
ethoxide = Chem.MolFromSmiles("[CH3CH2O-]")
sodium = Chem.MolFromSmiles("[Na+]")
bromobutane = Chem.MolFromSmiles("CCCCBr")

# Apply the reaction
products = reaction.RunReactants((ethoxide, sodium, bromobutane))

# Generate and save reaction image
img = Draw.ReactionToImage(reaction)
img.show()  # Display the image
img.save("williamson_ether_synthesis.png")  # Save the image
