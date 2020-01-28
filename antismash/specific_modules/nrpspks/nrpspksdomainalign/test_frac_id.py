import pyximport; pyximport.install()
import fraction_id_calc as fic

a = 'ABCDEFG----ABC--ABC'
b = 'GBCDEFGHI--ABWWWWBC'

print str(fic.ficalc(a, b))
