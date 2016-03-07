from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

# 1) directorio donde se encuentran los ficheros con coordenadas de moldes/templates,
# con extension .pdb,.atm,.ent
env.io.atom_files_directory = './template/'

# 2) prepara el modelado
a = automodel(env,
              alnfile  = './template/d1.ali',  # fichero con el alineamiento
              knowns   = '3TBG_A',              # nombre del template como aparece en alnfile
              sequence = 'd2j0da1 ',             # nombre de secuencia problema como aparece en alnfile
              assess_methods=(assess.DOPE))

a.starting_model= 1                           # define cuantos modelos diferentes quieres
a.ending_model  = 2

# 3) accion!
a.make()
