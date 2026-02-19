
# Python code that chooses  N proteins from a mfasta file
import random

random.seed(20250307)

def rinkinys(mfasta, N=100):
  f=open(mfasta, "r")
  eilutes=f.readlines()
  nuo = random.randint(0, len(eilutes)-N*30)
  f.close()
  g=open(mfasta+f'.{N}','w')
  prad=0

  while eilutes[prad][0]!='>':
     prad+=1
  N+=1
  atkarpa=eilutes[prad:nuo+N*30]
  for e in atkarpa:
    if e[0]=='>':
      N-=1
      if N==0:
        break
    g.write(e)
  g.close()
  return

rinkinys('viral.1.faa', N=100)
