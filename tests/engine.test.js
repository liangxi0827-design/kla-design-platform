// kLa 设计台计算引擎参考算例测试
// 运行：node tests/engine.test.js （需要 playwright；页面在本地以 file:// 打开，外部资源全部拦截）
const {chromium}=require('playwright');
const path=require('path');
(async()=>{
  const b=await chromium.launch(); const pg=await b.newPage();
  await pg.route('**/*',r=>r.request().url().startsWith('file:')?r.continue():r.abort());
  await pg.goto('file://'+path.resolve(__dirname,'../index.html')); await pg.waitForTimeout(800);
  const results=await pg.evaluate(()=>{
    const E=window.__klaEngine, out=[];
    const T=(name,ok,detail)=>out.push({name,ok:!!ok,detail});
    const base=Object.assign({},E.BASE,{model:"vr_c",fM:1,fA:1,fmed:1,tcorr:false,fsalt:1,fco2:1,fmix:1});
    const close=(a,b,rel)=>Math.abs(a-b)<=rel*Math.max(Math.abs(a),Math.abs(b),1e-12);

    // 1 氧溶解度：Benson–Krause，20 °C 纯水、1 atm 空气约 9.09 mg/L
    const c20=E.csWater(20); T("C*(20 °C) = 9.09 mg/L ±0.3%",close(c20,9.09,0.003),c20.toFixed(3)+" mg/L");
    const c37=E.csWater(37); T("C*(37 °C) 在 6.6–6.8 mg/L",c37>6.6&&c37<6.8,c37.toFixed(3)+" mg/L");

    // 2 几何体积闭合：由 T、H 反算体积 = 输入体积
    for(const head of ["ell","tori","flat"]){ const p=Object.assign({},base,{vol:3000,ht:1.5,head}); const g=E.geometry(p);
      const k={ell:[0.25,Math.PI/24],tori:[0.169,0.0809],flat:[0,0]}[head]; const V=Math.PI/4*g.T**2*(g.H-k[0]*g.T)+k[1]*g.T**3;
      T(`体积闭合（${head}）`,close(V,3,1e-9),`反算 ${(V*1000).toFixed(3)} L`); }

    // 3 功率：单叶轮 P0 = Np·ρ·N³·D⁵
    { const p=Object.assign({},base,{nImp:1,impB:"rushton",npB:5.5,rpm:60,microOn:false,macroOn:false}); const r=E.compute(p);
      const ref=5.5*p.rho*1**3*r.D**5; T("P0 = Np·ρ·N³·D⁵",close(r.P0,ref,1e-9),`${r.P0.toFixed(2)} vs ${ref.toFixed(2)} W`); }

    // 4 Van’t Riet 复现：单分布器、倍数 1、无校正时 kLa = 0.026·(Pg/V)^0.4·vs^0.5（vs 为工况体积）
    { const p=Object.assign({},base,{macroOn:false,microOn:true,vvmM:0.01,rpm:80}); const r=E.compute(p);
      const ref=0.026*Math.pow(r.PgV,0.4)*Math.pow(r.vs,0.5)*3600; T("Van’t Riet 复现",close(r.kla,ref,1e-9),`${r.kla.toFixed(3)} vs ${ref.toFixed(3)} h⁻¹`); }

    // 5 标况→工况：vs 按 (T/Tref)·(1 atm/P̄) 换算
    { const p=Object.assign({},base,{macroOn:false,vvmM:0.01}); const r=E.compute(p); const A=Math.PI*r.T**2/4;
      const ref=r.QM*((p.temp+273.15)/273.15)*(1.01325/r.Pm)/A; T("vs 用工况体积流量",close(r.vs,ref,1e-9),`${(r.vs*1000).toFixed(4)} mm/s`); }

    // 6 氧平衡：任何工况下 OTR·V ≤ 通入 O₂（利用率 ≤ 100%）
    { let worst=0; for(const vol of [200,2000,12500,20000]) for(const vv of [0.001,0.005,0.02]) for(const rpm of [30,80,150]){
        const r=E.compute(Object.assign({},base,{vol,vvmM:vv,vvmA:vv,rpm,xo2M:1}),{}); worst=Math.max(worst,r.o2util); }
      T("氧平衡：O₂ 利用率 ≤ 100%",worst<=1+1e-9,`最大利用率 ${(worst*100).toFixed(1)}%`); }

    // 7 极限一致性：通气量很大（NTU→0）时，计气相消耗的 OTR → kLa·(C*−C_L)
    { const p=Object.assign({},base,{macroOn:false,vvmM:0.5,rpm:40}); const r=E.compute(p);
      T("大气量极限：OTR → kLa(C*−C_L)",close(r.OTR,r.OTRnd,0.02),`${r.OTR.toFixed(3)} vs ${r.OTRnd.toFixed(3)} mmol/L/h`); }

    // 8 CO₂ 物料平衡：由 pCO₂ 反算的去除量 = CER·V
    { const p=Object.assign({},base,{vvmM:0.01,vvmA:0.02}); const r=E.compute(p); const Tref=273.15, G=q=>q*101325/(8.314*Tref);
      const He=34*Math.exp(2400*(1/(p.temp+273.15)-1/298.15)), PmA=r.Pm/1.01325, ks=r.kla/3600, V=r.V;
      const rem=[[G(r.QM),r.sM],[G(r.QA),r.sA]].reduce((u,[Gi,s])=>u+(Gi>0?Gi*(1-Math.exp(-(p.kco2*ks*s)*V*He*PmA/Gi)):0),0);
      const removed=(r.pco2/760)/PmA*rem; T("CO₂ 物料平衡闭合",close(removed,r.CER*V,1e-9),`${removed.toExponential(4)} vs ${(r.CER*V).toExponential(4)} mol/s`); }

    // 9 孔速手算：Q 实际（罐底压力、培养温度）/ 孔总面积
    { const p=Object.assign({},base,{macroOn:true,vvmA:0.02,holesA:60,holeD:1.5}); const r=E.compute(p);
      const q=r.QA*((p.temp+273.15)/273.15)*(1.01325/r.Pbot), ref=q/(60*Math.PI*0.0015**2/4); T("出口气速手算",close(r.gevA,ref,1e-9),`${r.gevA.toFixed(2)} m/s`); }

    // 10 OUR 单位：qO₂[pmol/细胞/h]·VCD[10⁶/mL] = mmol/L/h
    { const r=E.compute(Object.assign({},base,{qo2:0.3,vcd:20})); T("OUR = qO₂·VCD",close(r.OUR,6,1e-12),`${r.OUR} mmol/L/h`); }

    // 11 单调性：kLa 随转速、通气增加
    { const k=(rpm,v)=>E.compute(Object.assign({},base,{rpm,vvmM:v})).kla; T("kLa 随转速、通气单调增加",k(40,0.01)<k(80,0.01)&&k(80,0.01)<k(80,0.02),""); }

    // 12 判据三态：估算类判据在限值 ±20% 内给“证据不足”
    { const r0=E.compute(base), lim=r0.theta*1.1, r=E.compute(Object.assign({},base,{mixMax:lim})); const c=E.evaluate(r).find(x=>x.lbl.startsWith("混合时间"));
      T("三态判据：估算值距限值 <20% → 证据不足",c&&c.st==="warn",`θ=${r.theta.toFixed(0)} s，限值 ${lim.toFixed(0)} s → ${c&&c.st}`);
      const r2=E.compute(Object.assign({},base,{mixMax:r0.theta*3})); const c2=E.evaluate(r2).find(x=>x.lbl.startsWith("混合时间"));
      T("三态判据：远离限值 → 满足",c2&&c2.st==="ok",`→ ${c2&&c2.st}`); }
    return out;
  });
  let pass=0; results.forEach(t=>{ console.log(`${t.ok?"PASS":"FAIL"}  ${t.name}  ${t.detail||""}`); if(t.ok) pass++; });
  console.log(`\n${pass}/${results.length} 通过`);
  await b.close(); process.exit(pass===results.length?0:1);
})();
