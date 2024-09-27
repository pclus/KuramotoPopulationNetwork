init=run('oa2',ICP=[1,2,3],IPS=-2,NMX=100000,PAR={'K': 6.0, 'p' : 0.9, 'alpha' : 1.2})

init=run('oa2',ICP=[1,2,3],IPS=-2,NMX=100000,PAR={'K': 7.0, 'p' : 0.9, 'alpha' : 1.2})
ic=init(201)
init=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2])

hb=init('HB1')
c1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=50000,ISW=1,DSMAX=0.01,NTST=200,NCOL=7, STOP=['BP1'])

init=run('oa2',ICP=[1,2,3],IPS=-2,NMX=100000,PAR={'K': 10.0, 'p' : 0.9, 'alpha' : 1.2})
ic=init(201)
init=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],DS="-")

hb=init('HB1')
c1=run(hb,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=50000,ISW=1,DSMAX=0.01,NTST=200,NCOL=7, UZSTOP={'K' : 9.3})
unstable = c1('UZ1')
save(unstable,'unstableLC')
c1=run(unstable,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=50000,ISW=1,DSMAX=0.01,NTST=200,NCOL=7, UZSTOP={'K' : 9.3})

stable = c1('UZ1')
save(stable,'stableLC')


init=run('oa2',ICP=[1,2,3],IPS=-2,NMX=100000,PAR={'K': 7.0, 'p' : 0.9, 'alpha' : 1.2})
ic=init(201)
init=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'K':9.3})


