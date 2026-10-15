from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).parent
YELLOW='#FFD500'; DARK='#202020'; GRAY='#777777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def save(fig,name):
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/f'{name}.png',bbox_inches='tight',dpi=160,facecolor='white')
    plt.close(fig)

fig,ax=plt.subplots(figsize=(3,4))
ax.axis('off'); ax.set(xlim=(0,1),ylim=(0,1))
for i in range(100):
    ax.add_patch(plt.Rectangle((.32,.12+i*.007),.36,.0075,color=plt.cm.inferno(.25+i*.006)))
ax.plot([.32,.32],[.12,.82],c=DARK,lw=2); ax.plot([.68,.68],[.12,.82],c=DARK,lw=2)
ax.annotate('',(.5,.18),(.5,.77),arrowprops={'arrowstyle':'->','lw':2,'color':'white'})
ax.text(.5,.91,'He + vapour',ha='center',weight='bold')
ax.text(.5,.035,'Cooling column',ha='center')
ax.text(.76,.72,'hot',ha='left'); ax.text(.76,.22,'cold',ha='left')
save(fig,'title-column')

fig,ax=plt.subplots(figsize=(6,2.1)); ax.axis('off'); ax.set(xlim=(0,10),ylim=(0,3))
for x,txt,col in [(0.1,'Cell gas\npartial pressure', '#eeeeee'),(3.7,'Diffusion\nresistance',YELLOW),(7.2,'Wall\ninventory','#eeeeee')]:
    ax.add_patch(plt.Rectangle((x,.6),2.6,1.7,facecolor=col,edgecolor=DARK))
    ax.text(x+1.3,1.45,txt,ha='center',va='center')
for a,b in [(2.8,3.6),(6.4,7.1)]: ax.annotate('',(b,1.4),(a,1.4),arrowprops={'arrowstyle':'<->','lw':2})
save(fig,'wall-law')

fig,ax=plt.subplots(figsize=(6,3)); ax.axis('off'); ax.set(xlim=(0,10),ylim=(0,5))
for x,y,w,h,txt,col in [(0,2.7,10,1.5,'Pb   PbI   PbI₂   Bi   Bi₂   BiI   BiI₃   I   I₂',YELLOW),(0,.1,4.3,1.5,'Gas equilibrium\nPb / Bi / I preserved','#eeeeee'),(5.7,.1,4.3,1.5,'Shared wall stores\nPbI₂   BiI₃   Pb   Bi','#eeeeee')]:
    ax.add_patch(plt.Rectangle((x,y),w,h,facecolor=col,edgecolor=DARK)); ax.text(x+w/2,y+h/2,txt,ha='center',va='center')
for x in [2.2,7.8]: ax.annotate('',(x,1.7),(x,2.6),arrowprops={'arrowstyle':'<->','lw':2})
save(fig,'chemistry')

fig,ax=plt.subplots(figsize=(6.1,3.8))
x=np.arange(3)
ax.plot(x,[15.35,8.95,1.98],'-o',color=DARK,lw=2,ms=8,label='limitedLinear')
ax.plot(x,[10.47,10.58,2.13],'--s',color='#aa8500',lw=2,ms=7,label='Upwind')
ax.set_xticks(x,['8,192 → 18,432','18,432 → 41,472','41,472 → 93,312'],fontsize=9)
ax.set(ylabel='Iodine profile change (%)',xlabel='Successive mesh pairs (cells)',ylim=(0,18))
ax.legend(frameon=False); ax.grid(axis='y',alpha=.2)
ax.text(2,3.1,'1.98%',ha='center',weight='bold')
save(fig,'numerics')

fig,ax=plt.subplots(figsize=(6.1,3.8))
x=np.arange(5); w=.23
for offset,values,color,label in [(-w,[26,19,16,-11,-14],YELLOW,'This work · RMS 18 K'),(0,[-39,-10,-1,-60,-30],DARK,'Published temperature · 35 K'),(w,[27,37,34,-50,-15],'#aaaaaa','Published vapour pressure · 35 K')]:
    ax.bar(x+offset,values,w,color=color,label=label,edgecolor='white',linewidth=.5)
ax.axhspan(-15,15,color='#dddddd',alpha=.25,zorder=0); ax.axhline(0,color=DARK,lw=.7)
ax.set_xticks(x,['106 / 6','106 / 60','106 / 300','5 / 60','540 / 60'])
ax.set(ylabel='Onset minus measurement (K)',xlabel='Flow (mL/min) / exposure (s)',ylim=(-68,63))
ax.legend(fontsize=8.5,frameon=False,loc='upper left'); save(fig,'onset-errors')

fig,ax=plt.subplots(figsize=(5.9,3.7))
values=[94.58,2.01]
bars=ax.bar([0,1],values,color=[YELLOW,DARK],width=.52)
ax.set_xticks([0,1],['Baseline BiI data','BiI +29 kJ/mol']); ax.set(ylabel='BiI₃ share of deposited iodide iodine (%)',ylim=(0,108))
for b,v in zip(bars,values): ax.text(b.get_x()+b.get_width()/2,v+3,f'{v:.2f}%',ha='center',weight='bold',fontsize=15)
save(fig,'partition')

fig,ax=plt.subplots(figsize=(6,3.7))
values=[6.13,77.80,73.26]
bars=ax.barh([2,1,0],values,color=[DARK,YELLOW,'#c5a700'],height=.55)
ax.set_yticks([2,1,0],['Original input','−4.10 cm shift','−4.65 cm shift']); ax.set(xlim=(0,100),xlabel='Normalised scan overlap (%)')
for b,v in zip(bars,values): ax.text(v+2,b.get_y()+b.get_height()/2,f'{v:.2f}%',va='center',weight='bold')
save(fig,'silica')

