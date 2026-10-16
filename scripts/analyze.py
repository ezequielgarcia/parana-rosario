import json, datetime as dt, statistics as st

def load(i):
    d=json.load(open('/tmp/claude-1000/-home-zeta-parana/5456bbce-43b9-4fdc-97d9-fd46dd19c387/scratchpad/dl/ser_%d.json'%i))
    s={}
    for o in d['observaciones']:
        v=o.get('valor')
        if v is None: continue
        day=dt.date.fromisoformat(o['timestart'][:10])
        s[day]=float(v)
    return s

C=load(19); R=load(34); SF=load(30)
start=min(min(C),min(R)); end=max(max(C),max(R))
days=[start+dt.timedelta(d) for d in range((end-start).days+1)]

def series(s): return [s.get(d) for d in days]
c=series(C); r=series(R)

# 1) cross-correlation of daily increments
def xcorr(a,b,lag):
    xs=[];ys=[]
    for i in range(1,len(a)-lag):
        if a[i] is None or a[i-1] is None: continue
        j=i+lag
        if j<1 or j>=len(b) or b[j] is None or b[j-1] is None: continue
        xs.append(a[i]-a[i-1]); ys.append(b[j]-b[j-1])
    n=len(xs)
    mx=sum(xs)/n; my=sum(ys)/n
    num=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    den=(sum((x-mx)**2 for x in xs)*sum((y-my)**2 for y in ys))**0.5
    return num/den, n

print("lag  corr(dCorrientes, dRosario)   n")
best=None
for lag in range(0,25):
    cc,n=xcorr(c,r,lag)
    print(f"{lag:3d}  {cc:+.3f}  {n}")
    if best is None or cc>best[1]: best=(lag,cc)
print("best lag (increments):",best)

# 2) event matching: Corrientes peaks above 6.0 m, separated by >=40 days
peaks=[]
for i in range(10,len(c)-30):
    if c[i] is None or c[i]<6.0: continue
    win=[x for x in c[max(0,i-12):i+13] if x is not None]
    if c[i]==max(win) and all((c[i]>=x) for x in win):
        peaks.append(i)
# dedupe
ded=[]
for i in peaks:
    if ded and i-ded[-1]<40:
        if c[i]>c[ded[-1]]: ded[-1]=i
    else: ded.append(i)

print("\nCorrientes peak -> Rosario response")
print("date_C   C_peak   date_R   R_peak   lag_d   dC_rise30  dR_rise")
rows=[]
for i in ded:
    # Rosario peak within i-5 .. i+30
    seg=[(j,r[j]) for j in range(i-5,min(i+31,len(r))) if r[j] is not None]
    if not seg: continue
    j,rp=max(seg,key=lambda t:t[1])
    # require the Rosario max not be at the window edge
    base_c=[x for x in c[max(0,i-40):i-25] if x is not None]
    base_r=[x for x in r[max(0,i-40):i-25] if x is not None]
    if not base_c or not base_r: continue
    bc=st.median(base_c); br=st.median(base_r)
    rows.append((days[i],c[i],days[j],rp,(days[j]-days[i]).days,c[i]-bc,rp-br))
    print(f"{days[i]} {c[i]:6.2f}  {days[j]} {rp:6.2f} {(days[j]-days[i]).days:5d}  {c[i]-bc:+8.2f} {rp-br:+8.2f}")

# regression R_peak = a + b*C_peak
xs=[x[1] for x in rows]; ys=[x[3] for x in rows]
n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
b=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
a=my-b*mx
res=[y-(a+b*x) for x,y in zip(xs,ys)]
sd=(sum(e*e for e in res)/(n-2))**0.5
print(f"\nn={n}  R_peak = {a:.2f} + {b:.3f}*C_peak   resid sd={sd:.2f} m")
lags=[x[4] for x in rows]
print("lags: median %.0f mean %.1f min %d max %d" % (st.median(lags), sum(lags)/len(lags), min(lags), max(lags)))

for cp in (6.75,7.31,7.55,8.0):
    print(f"  C_peak {cp:.2f} -> R_peak {a+b*cp:.2f}  (+/-{sd:.2f})")

# 3) rise-based: Rosario rise vs Corrientes rise
xs2=[x[5] for x in rows]; ys2=[x[6] for x in rows]
mx2=sum(xs2)/n; my2=sum(ys2)/n
b2=sum((x-mx2)*(y-my2) for x,y in zip(xs2,ys2))/sum((x-mx2)**2 for x in xs2)
a2=my2-b2*mx2
res2=[y-(a2+b2*x) for x,y in zip(xs2,ys2)]
sd2=(sum(e*e for e in res2)/(n-2))**0.5
print(f"\nR_rise = {a2:.2f} + {b2:.3f}*C_rise  resid sd={sd2:.2f}")
