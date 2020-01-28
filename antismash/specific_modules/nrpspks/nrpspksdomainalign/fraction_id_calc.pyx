def ficalc(a, b):
 m, t = 0, 0
 lista, listb = list(a), list(b)
 for i in range(0,len(a)):
  if lista[i] == '-':
   continue
  if listb[i] == '-':
   continue
  t += 1
  if lista[i] == listb[i]:
   m += 1
 if t == 0:
  return 0
 else:
  return round(float(m)/t, 4)