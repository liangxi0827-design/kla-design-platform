// 批次动态模拟测试：node tests/sim.test.js
const {chromium}=require('playwright'); const path=require('path');
(async()=>{ const b=await chromium.launch(); const pg=await b.newPage();
  await pg.route('**/*',r=>r.request().url().startsWith('file:')?r.continue():r.abort());
  await pg.goto('file://'+path.resolve(__dirname,'../index.html')); await pg.waitForTimeout(800);
  const res=await pg.evaluate(()=>{ const E=window.__klaSim, out=[], T=(n,ok,d)=>out.push({n,ok:!!ok,d});
    const base=Object.assign({},E.BASE,E.PRESETS.xdr2000.v), P0=Object.assign({},E.SIMP);
    // 1 无生长、无补料：活细胞按 kd 指数衰减
    { const P=Object.assign({},P0,{muMax:0,feedPct:0,kdA:0,kdS:0,feedStart:99,days:2}); const r=E.simulate(base,P,{kineticOnly:true});
      const exp=P.X0*Math.exp(-P.kd0*48), got=r.final.Xv; T("μ=0 时 Xv 按 kd 衰减",Math.abs(got/exp-1)<0.01,`${got.toFixed(4)} vs ${exp.toFixed(4)}`); }
    // 2 活细胞 + 死细胞总量守恒（无补料稀释时）
    { const P=Object.assign({},P0,{feedPct:0,feedStart:99,days:6}); const r=E.simulate(base,P,{kineticOnly:true});
      const ok=r.out.every((o,i)=>i===0||o.Xv+o.Xd>=r.out[i-1].Xv+r.out[i-1].Xd-1e-9); T("无补料时总细胞数不减少",ok,""); }
    // 3 体积不超过工作体积
    { const r=E.simulate(base,Object.assign({},P0,{feedPct:10}),{kineticOnly:true}); const mx=Math.max(...r.out.map(o=>o.V)); T("体积 ≤ 工作体积",mx<=base.vol+1e-6,`最大 ${mx.toFixed(1)} L`); }
    // 4 供氧充足时 DO 维持设定
    { const r=E.simulate(base,Object.assign({},P0,{days:6})); const mn=Math.min(...r.out.map(o=>o.DO)); T("供氧充足时 DO = 设定",Math.abs(mn-base.doset)<0.5,`最低 ${mn.toFixed(1)}%`); }
    // 5 只通空气且不允许调整时，高密度下 DO 下降
    { E.SIMP.useVM=false; E.SIMP.useRpm=false; const xm=E.SIMP.xMax; E.SIMP.xMax=0.21;
      const r=E.simulate(Object.assign({},base,{xo2M:0.21}),Object.assign({},P0,{days:12}));
      E.SIMP.useVM=true; E.SIMP.useRpm=true; E.SIMP.xMax=xm;
      const mn=Math.min(...r.out.map(o=>o.DO)); T("控制受限时 DO 跌破设定",mn<base.doset-5,`最低 ${mn.toFixed(1)}%`); }
    // 6 合成数据参数回收（无噪声）
    { const Ptrue=Object.assign({},P0,{muMax:P0.muMax*0.85,qP:P0.qP*1.2}); const sim=E.simulate(base,Ptrue,{kineticOnly:true});
      const data=[1,2,3,4,5,6,7,8,9,10,11,12,13,14].map(d=>{ const o=E.simAt(sim,d); return {day:d,Xv:o.Xv,via:o.via,G:o.G,L:o.L,P:o.P}; });
      const keys=["muMax","qP"]; const f=z=>{ const P=Object.assign({},P0); keys.forEach((k,i)=>P[k]=P0[k]*Math.exp(z[i])); return E.batchObj(E.simulate(base,P,{kineticOnly:true}),data); };
      const r=E.nelderMead(f,[0,0],150); const mu=P0.muMax*Math.exp(r.x[0]), qp=P0.qP*Math.exp(r.x[1]);
      T("参数回收：μmax、qP 误差 < 5%",Math.abs(mu/Ptrue.muMax-1)<0.05&&Math.abs(qp/Ptrue.qP-1)<0.05,`μmax ${mu.toFixed(4)}/${Ptrue.muMax.toFixed(4)}，qP ${qp.toFixed(2)}/${Ptrue.qP.toFixed(2)}`); }
    // 7 OUR = qO₂·Xv（降温后乘系数）
    { const r=E.simulate(base,P0,{kineticOnly:true}); const o=E.simAt(r,10); const exp=base.qo2*P0.qo2Shift*o.Xv; T("OUR = qO₂·VCD",Math.abs(o.OUR/exp-1)<1e-9,`${o.OUR.toFixed(3)}`); }
    return out; });
  let p=0; res.forEach(t=>{ console.log(`${t.ok?"PASS":"FAIL"}  ${t.n}  ${t.d}`); if(t.ok) p++; }); console.log(`\n${p}/${res.length} 通过`); await b.close(); process.exit(p===res.length?0:1); })();
