// CFD 降阶模型测试：从 index.html 提取求解核心，验证守恒、对称性与关键物理量的量级
const fs=require('fs'), path=require('path');
const html=fs.readFileSync(path.join(__dirname,'..','index.html'),'utf8');
const a=html.indexOf('function CFDCORE(){'), b=html.indexOf('return {run,bubbleUt,makeGrid,lbmRun};\n}',a);
const CFDCORE=new Function(html.slice(a,b)+'return {run,bubbleUt,makeGrid,lbmRun};\n}\nreturn CFDCORE;')();
const C=CFDCORE();
let pass=0, fail=0; const ok=(name,cond,info)=>{ if(cond){pass++; console.log('✓',name,info||'');} else {fail++; console.log('✗',name,info||'');} };

// 1. 轴对称 LBM：Hagen–Poiseuille 管流中心速度 = gR²/(4ν)
{ const nx=24, cfg={T:1,Hl:4/24,hd:0,imps:[],plumes:[],Pimp:0,Pgas:0};
  for(const axi of [true,false]){ const g=C.makeGrid(cfg,nx); g.periodic=true; const nu=0.05, gF=1e-6;
    const L=C.lbmRun(cfg,g,{steps:12000,avgFrom:11000,nu0:nu,Cs:0,axi,bodyY:gF},null); const R=nx/2;
    let err=0,mx=0; for(let i=0;i<nx;i++){ const x=i+0.5-R, ex=gF*(R*R-x*x)/((axi?4:2)*nu); err=Math.max(err,Math.abs(L.uy[nx+i]-ex)); mx=Math.max(mx,ex); }
    ok(`${axi?'轴对称':'平面'} Poiseuille 解析解误差 < 3%`,err/mx<0.03,(err/mx*100).toFixed(1)+'%'); } }

// 2–7. 3000 L 不锈钢双斜叶算例
const T=1.392,H=2.088,D=0.557,R=T/2;
const cfg={T,Hl:H,hd:0.25*T,tilt:0,N:1,rho:1005,nu:0.9e-6,V:3,Pimp:140,Pgas:7.5,sym:true,axi:false,plumeForce:true,
 imps:[{x:0,y:0.418,D,Nq:0.79,angle:35,up:false,radial:false,type:'pbt4'},{x:0,y:1.253,D,Nq:0.79,angle:35,up:false,radial:false,type:'pbt4'}],
 plumes:[{x:0.234,y:0.25,w:0.03,q:1},{x:-0.234,y:0.25,w:0.03,q:1}],
 gas:[{name:'微泡',Q:0.00028,d0:0.0009,yEff:0.45,src:[{x:0.25,y:0.25,w:0.05},{x:-0.25,y:0.25,w:0.05}]},{name:'大泡',Q:0.0009,d0:0.004,yEff:0.2,src:[{x:0.234,y:0.26,w:0.02},{x:-0.234,y:0.26,w:0.02}]}],
 sigma:0.06,DO2:2.2e-9,fmed:1,HO2:6.7/32/0.2095,Phead:1.31325,OUR:0.29*30/3600,cL:0.4*0.215,cAirHead:0.215,
 inject:{x:0.5*R,y:H*0.97},probes:[{x:0.6*R,y:0.12*H},{x:-0.6*R,y:0.5*H},{x:0.6*R,y:0.9*H}],tMixMax:300,tDO:80,doProbe:{x:-0.85*R,y:0.25*H},doLow:20,dCell:18e-6,nx:36,fineCirc:1.6};
const res=C.run(cfg), s=res.stats, gr=res.grid, F=res.fields;
let e=0; for(let c=0;c<gr.wv.length;c++) if(!gr.solid[c]) e+=gr.wv[c]*F.eps[c]*cfg.rho;
ok('能量守恒：∫ε·ρ dV = Pg + P_gas（±0.5%）',Math.abs(e/(cfg.Pimp+cfg.Pgas)-1)<0.005,(e).toFixed(2)+' W');
let asym=0, mx=0; for(let j=0;j<gr.ny;j++) for(let i=0;i<gr.nx/2;i++){ const p=j*gr.nx+i, q=j*gr.nx+gr.nx-1-i; asym=Math.max(asym,Math.abs(F.ux[p]+F.ux[q]),Math.abs(F.uy[p]-F.uy[q])); mx=Math.max(mx,F.spd[p]); }
ok('中心安装平均流场关于轴线对称',asym<1e-6*mx+1e-9);
ok('预测泵送数在厂商 Nq 的 ±50% 内',Math.abs(s.NqPred/0.79-1)<0.5,s.NqPred.toFixed(2)+' vs 0.79');
ok('平均流速为叶尖速度的 5–45%',s.uMean/s.Utip>0.05&&s.uMean/s.Utip<0.45,(s.uMean/s.Utip*100).toFixed(0)+'%');
const theta=5.9*T**(2/3)*Math.pow(s.epsMean,-1/3)*Math.pow(0.4,-1/3)*Math.pow(1.5,2.43);
ok('CFD 混合时间与 Grenville 关联式同一量级（0.33–3×）',s.t95/theta>0.33&&s.t95/theta<3,s.t95.toFixed(0)+' s vs '+theta.toFixed(0)+' s');
ok('混合终态：最后一帧示踪剂 CoV < 5%',res.mix.cov[res.mix.cov.length-1][1]<0.05);
const pr=res.grid; // DO 探头闭环
{ const i=Math.floor((cfg.doProbe.x+R)/pr.dx), j=Math.floor(cfg.doProbe.y/pr.dx); const v=F.DO[j*pr.nx+i];
  ok('DO 控制：探头处 DO = 设定值 ±2%',Math.abs(v-40)<2,v.toFixed(1)+'%'); }
ok('DO 非负且有限',s.doMin>=0&&isFinite(s.doMax));
ok('平均气含率 ≈ Σ Q·τ / V（拉格朗日一致性 ±15%）',(()=>{ const est=cfg.gas.reduce((a,cl,k)=>a+cl.Q*s.gas[k].tRes,0)/s.epsMean/0+0; const hv=cfg.gas.reduce((a,cl,k)=>a+cl.Q*s.gas[k].tRes,0); let vol=0; for(let c=0;c<gr.wv.length;c++) if(!gr.solid[c]) vol+=gr.wv[c]; return Math.abs(s.holdMean/(hv/vol)-1)<0.15; })());
ok('Kolmogorov 尺度量级 10–1000 µm',s.etaMin>1e-5&&s.etaMin<1e-3,(s.etaMin*1e6).toFixed(0)+' µm');
// 8. 气泡终端速度：1 mm ≈ 10 cm/s，4 mm ≈ 22 cm/s（污染水中刚性气泡）
ok('气泡终端速度 1 mm',Math.abs(C.bubbleUt(1e-3,1e-6,1000,0.06)-0.105)<0.02,(C.bubbleUt(1e-3,1e-6,1000,0.06)*100).toFixed(1)+' cm/s');
ok('气泡终端速度 4 mm',Math.abs(C.bubbleUt(4e-3,1e-6,1000,0.06)-0.22)<0.03,(C.bubbleUt(4e-3,1e-6,1000,0.06)*100).toFixed(1)+' cm/s');
console.log(`\n${pass}/${pass+fail} 通过`); process.exit(fail?1:0);
