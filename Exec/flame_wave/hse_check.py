#Declare Constants
g=2.0e14

#read output.dat for network data used in Prob_2d.f90
prob=open('output.dat')
prob_r=[]
prob_dpdr=[]
prob_rhog=[]
n=0
for i in prob:
  columns=i.split(	)
  if len(columns)==4:
    prob_r.append(float(columns[0]))
    prob_dpdr.append(float(columns[1]))
    prob_rhog.append(float(columns[2]))
  if prob_dpdr[n]==0:
    break
  else:
    n=n+1


#read data from flame_wave.hse
net=open('flame_wave.hse')
net_r=[]
net_p=[]
net_dpdr=[]
net_rhog=[]
n=0
for i in net:
  columns=i.split(          )
  if len(columns)==8:
    net_r.append(float(columns[0]))
    net_p.append(float(columns[3]))
    net_rhog.append(float(columns[1])*g)

#crunch the numbers
n=0
while n<=len(prob_r):
  if n==0:
    net_dpdr.append(0)
  else:
    net_dpdr.append((abs(net_p[n]-net_p[n-1]))/(net_r[n]-net_r[n-1]))
  n=n+1

pe_dpdr=[]
pe_rhog=[]
pe_net=[]
n=0
while n+1<len(prob_r):
  pe_dpdr.append(abs(net_dpdr[n+1]-prob_dpdr[n])/net_dpdr[n+1])
  pe_rhog.append(abs(net_rhog[n+1]-prob_rhog[n])/net_rhog[n+1])
  pe_net.append(abs(net_dpdr[n+1]-net_rhog[n+1])/net_dpdr[n+1])
  n=n+1

hse=open('hse.dat','w+')
n=0
while n+1<len(prob_r):
  hse.write(str(prob_r[n]))
  hse.write('	')
  hse.write(str(prob_dpdr[n]))
  hse.write('	')
  hse.write(str(net_dpdr[n+1]))
  hse.write('	')
  hse.write(str(pe_dpdr[n]))
  hse.write('	')
  hse.write(str(prob_rhog[n]))
  hse.write('	')
  hse.write(str(net_rhog[n+1]))
  hse.write('	')
  hse.write(str(pe_rhog[n]))
  hse.write('	')
  hse.write(str(pe_net[n]))
  hse.write('\n')
  n=n+1
