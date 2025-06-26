function EIoverL3(E,I,L){
    const factor=E*I/Math.pow(L,3);
    return [
        [12*factor, 6*L*factor, -12*factor, 6*L*factor],
        [6*L*factor, 4*L*L*factor, -6*L*factor, 2*L*L*factor],
        [-12*factor, -6*L*factor, 12*factor, -6*L*factor],
        [6*L*factor, 2*L*L*factor, -6*L*factor, 4*L*L*factor]
    ];
}

function applyBC(K,F,fixed){
    const n=K.length;
    const free=[];
    for(let i=0;i<n;i++) if(!fixed.includes(i)) free.push(i);
    const Kmod=[],Fmod=[];
    for(let i=0;i<free.length;i++){
        Kmod[i]=[];
        for(let j=0;j<free.length;j++){
            Kmod[i][j]=K[free[i]][free[j]];
        }
        Fmod[i]=F[free[i]];
    }
    return {Kmod,Fmod,indices:free};
}

function gaussSolve(A,b){
    const n=b.length;
    for(let i=0;i<n;i++){
        let maxRow=i;
        for(let k=i+1;k<n;k++) if(Math.abs(A[k][i])>Math.abs(A[maxRow][i])) maxRow=k;
        [A[i],A[maxRow]]=[A[maxRow],A[i]];
        [b[i],b[maxRow]]=[b[maxRow],b[i]];
        const pivot=A[i][i];
        for(let j=i;j<n;j++) A[i][j]/=pivot;
        b[i]/=pivot;
        for(let k=i+1;k<n;k++){
            const factor=A[k][i];
            for(let j=i;j<n;j++) A[k][j]-=factor*A[i][j];
            b[k]-=factor*b[i];
        }
    }
    const x=new Array(n).fill(0);
    for(let i=n-1;i>=0;i--){
        let sum=b[i];
        for(let j=i+1;j<n;j++) sum-=A[i][j]*x[j];
        x[i]=sum;
    }
    return x;
}

function multiplyMatrixVector(A,x){
    const n=A.length; const m=x.length; const r=new Array(n).fill(0);
    for(let i=0;i<n;i++) for(let j=0;j<m;j++) r[i]+=A[i][j]*x[j];
    return r;
}

let crossSectionMap = {};
if (typeof require !== 'undefined' && typeof window === 'undefined') {
    try {
        const data = require('./steel_cross_sections.json');
        crossSectionMap = data.reduce((a,c)=>{a[c.profile]=c; return a;},{});
    } catch(e) {}
}
if (typeof window !== 'undefined' && window.steelCrossSectionsData) {
    crossSectionMap = window.steelCrossSectionsData;
}

function setCrossSections(data){
    crossSectionMap = data || {};
}

function getCrossSection(name){
    if((!crossSectionMap || Object.keys(crossSectionMap).length===0) &&
       typeof window !== 'undefined' && window.steelCrossSectionsData){
        crossSectionMap = window.steelCrossSectionsData;
    }
    return crossSectionMap ? crossSectionMap[name] : undefined;
}

function getSelfWeightLineLoads(spans, name){
    const cs = getCrossSection(name);
    if(!cs) return [];
    const w = cs.gk_kg_per_m * 9.81;
    const loads=[]; let pos=0;
    spans.forEach(s=>{
        const len=typeof s==='number'?s:s.length;
        loads.push({w,start:pos,end:pos+len});
        pos+=len;
    });
    return loads;
}

function computeInertia(cs){
    if(!cs) return 0;
    if(cs.series === 'sawn' || cs.series === 'glulam'){
        if(cs.Iy_m4 !== undefined) return cs.Iy_m4;
        const b = cs.b_mm/1000;
        const h = cs.h_mm/1000;
        return b*Math.pow(h,3)/12;
    }
    const h = cs.h_mm;
    const b = cs.b_mm;
    const tw = cs.tw_mm;
    const tf = cs.tf_mm;
    if(h===undefined || b===undefined || tw===undefined || tf===undefined)
        return cs.I_y !== undefined ? cs.I_y : (cs.Iz_m4 !== undefined ? cs.Iz_m4 : (cs.Iy_m4 !== undefined ? cs.Iy_m4 : 0));
    const hw = h - 2*tf;
    const Iweb = tw*Math.pow(hw,3)/12;
    const If = b*Math.pow(tf,3)/12;
    const d = hw/2 + tf/2;
    const Icalc = (2*(If + b*tf*d*d) + Iweb)/1e12; // to m^4
    return Icalc;
}

function computeWeakAxisInertia(cs){
    if(!cs) return 0;
    if(cs.Iz_m4 !== undefined) return cs.Iz_m4;
    if(cs.series === 'sawn' || cs.series === 'glulam'){
        const b = cs.b_mm/1000;
        const h = cs.h_mm/1000;
        return h*Math.pow(b,3)/12;
    }
    const h = cs.h_mm;
    const b = cs.b_mm;
    const tw = cs.tw_mm;
    const tf = cs.tf_mm;
    if(h===undefined || b===undefined || tw===undefined || tf===undefined)
        return 0;
    const hw = h - 2*tf;
    const Iweb = hw*Math.pow(tw,3)/12;
    const If = b*Math.pow(tf,3)/12;
    return (2*If + Iweb)/1e12;
}

function computeSectionDesign(name, opts){
    const cs = typeof name === 'string' ? getCrossSection(name) : name;
    if(!cs) return null;
    opts = opts || {};
    const material = opts.material || 'steel';
    const h = cs.h_mm/1000;
    const tw = cs.tw_mm/1000;
    const tf = cs.tf_mm/1000;
    let I = computeInertia(cs);
    const W = I/(h/2);

    if(material === 'timber'){
        const E = opts.E !== undefined ? opts.E : 11e9;
        const fm = opts.fm_k !== undefined ? opts.fm_k : 24e6;
        const fv = opts.fv_k !== undefined ? opts.fv_k : 4e6;
        const gammaM = opts.gammaM !== undefined ? opts.gammaM : 1.25;
        const EI = E*I;
        const Av = (cs.b_mm/1000)*h;
        const MRd = fm*W/gammaM;
        const VRd = fv*Av/gammaM;
        let MRdLBA = MRd;
        if(typeof opts.unbracedLength === 'number' && opts.unbracedLength > 0){
            const Lb = opts.unbracedLength;
            const b = cs.b_mm/1000;
            const hsec = cs.h_mm/1000;
            const fmd = fm/gammaM;
            const sigmaCrit = 0.78*b*b*E/(hsec*Lb);
            const lambdaRel = Math.sqrt(fmd/sigmaCrit);
            const kcrit = 1/(lambdaRel + Math.sqrt(lambdaRel*lambdaRel + 0.25));
            MRdLBA = kcrit*W*fmd;
        }
        return {EI, MRd, MRdLBA, VRd, W, gamma: gammaM, material: 'timber'};
    }

    const E = opts.E !== undefined ? opts.E : 210e9;
    const fy = opts.fy !== undefined ? opts.fy : 355e6;
    const gammaM0 = opts.gammaM0 !== undefined ? opts.gammaM0 : 1.0;
    const EI = E*I;
    const MRd = fy*W/gammaM0;
    const hw = h - 2*tf;
    const Av = hw*tw;
    const VRd = Av*fy/(Math.sqrt(3)*gammaM0);
    let MRdLBA = MRd;
    const b = cs.b_mm/1000;
    let It = ((2*b*Math.pow(tf,3))/3 + (hw*Math.pow(tw,3))/3);
    const y = h/2 - tf/2;
    let Iw = 2*(b*Math.pow(tf,3)/12)*Math.pow(y,2);
    let Mcr, chiLT;
    if(typeof opts.unbracedLength === 'number' && opts.unbracedLength > 0){
        const Lb = opts.unbracedLength;
        const Iz = computeWeakAxisInertia(cs);
        const G = opts.G !== undefined ? opts.G : 81e9;
        const C1 = opts.C1 !== undefined ? opts.C1 : 1.0;
        const C2 = opts.C2 !== undefined ? opts.C2 : 0.0;
        const C3 = opts.C3 !== undefined ? opts.C3 : 0.0;
        const kw = opts.kw !== undefined ? opts.kw : 1.0;
        const base = C1*Math.PI/Lb*Math.sqrt(E*Iz*G*It);
        const term2 = C2*(Math.PI*Math.PI*E*Iw)/(C1*C1*G*It*Lb*Lb);
        const term3 = C3*(Math.pow(Math.PI,4)*E*E*Iz*Iw)/(C1*C1*G*G*It*It*Math.pow(Lb,4));
        Mcr = base*Math.sqrt(kw + term2 + term3);
        const lambdaRel = Math.sqrt((fy*W/gammaM0)/Mcr);
        const alpha = 0.34;
        const phi = 0.5*(1 + alpha*(lambdaRel-0.2) + lambdaRel*lambdaRel);
        chiLT = 1/(phi + Math.sqrt(phi*phi + lambdaRel*lambdaRel));
        MRdLBA = chiLT*(fy*W/gammaM0);
    }
    return {EI, MRd, MRdLBA, VRd, W, gamma: gammaM0, material: 'steel', Iw, It, Mcr, chiLT};
}

function computeResults(state){
    if(!state.spans || state.spans.length===0) return null;
    const maxDist = state.maxNodeDist || 1;
    const E = state.E || 210e9;
    let I = state.I;
    if(I===undefined && state.section){
        const cs=getCrossSection(state.section);
        if(cs) I=computeInertia(cs);
    }
    if(I===undefined) I=1e-6;

    // Build nodes and supports
    const spanLengths = state.spans.map(s=>s.length);
    const boundaries=[0];
    let acc=0;
    spanLengths.forEach(len=>{acc+=len; boundaries.push(acc);});

    const posSet=new Set(boundaries);
    (state.pointLoads||[]).forEach(p=>posSet.add(p.x));
    (state.lineLoads||[]).forEach(l=>{posSet.add(l.start); posSet.add(l.end);});
    (state.extraNodePositions||[]).forEach(x=>posSet.add(x));

    const critical=Array.from(posSet).sort((a,b)=>a-b);
    const nodes=[critical[0]];
    for(let i=0;i<critical.length-1;i++){
        const start=critical[i];
        const end=critical[i+1];
        const segs=Math.ceil((end-start)/maxDist);
        for(let j=1;j<=segs;j++) nodes.push(start+(end-start)*j/segs);
    }

    const indexMap={};
    nodes.forEach((v,i)=>{indexMap[v.toFixed(6)]=i;});
    const supports=[];
    const leftSupport = state.leftSupport !== false;
    const rightSupport = state.rightSupport !== false;
    boundaries.forEach((b,i)=>{
        if(i===0 && !leftSupport) return;
        if(i===boundaries.length-1 && !rightSupport) return;
        const idx=indexMap[b.toFixed(6)];
        if(idx!==undefined) supports.push(idx);
    });
    const N=nodes.length;

    // global matrices
    const size=2*N;
    const K=[];
    const F=new Array(size).fill(0);
    for(let i=0;i<size;i++){K[i]=new Array(size).fill(0);}

    for(let e=0;e<N-1;e++){
        const x1=nodes[e];
        const x2=nodes[e+1];
        const L=x2-x1;
        const k=EIoverL3(E,I,L);
        const dofs=[2*e,2*e+1,2*(e+1),2*(e+1)+1];
        for(let i=0;i<4;i++){
            for(let j=0;j<4;j++){
                K[dofs[i]][dofs[j]]+=k[i][j];
            }
        }
    }

    // line loads
    (state.lineLoads||[]).forEach(l=>{
        for(let e=0;e<N-1;e++){
            const x1=nodes[e];
            const x2=nodes[e+1];
            const L=x2-x1;
            const overlap=Math.max(0,Math.min(l.end,x2)-Math.max(l.start,x1));
            if(overlap>0){
                const frac=overlap/L;
                const w=l.w;
                const f=[w*L*frac/2,w*L*L*frac/12,w*L*frac/2,-w*L*L*frac/12];
                const dofs=[2*e,2*e+1,2*(e+1),2*(e+1)+1];
                for(let i=0;i<4;i++) F[dofs[i]]+=f[i];
            }
        }
    });
    // point loads
    (state.pointLoads||[]).forEach(p=>{
        for(let e=0;e<N-1;e++){
            const x1=nodes[e];
            const x2=nodes[e+1];
            if(p.x>=x1 && p.x<=x2){
                const L=x2-x1;
                const a=p.x-x1, b=L-a;
                const P=p.P;
                const f=[P*b*b*(3*a+b)/Math.pow(L,3), P*a*b*b/Math.pow(L,2),
                         P*a*a*(3*b+a)/Math.pow(L,3), -P*a*a*b/Math.pow(L,2)];
                const dofs=[2*e,2*e+1,2*(e+1),2*(e+1)+1];
                for(let i=0;i<4;i++) F[dofs[i]]+=f[i];
                break;
            }
        }
    });

    // boundary conditions
    const fixed=[];
    supports.forEach(idx=>fixed.push(2*idx));

    const {Kmod,Fmod}=applyBC(K,F,fixed);
    const U=gaussSolve(Kmod,Fmod);
    const fullU=new Array(size).fill(0);
    let c=0;
    for(let i=0;i<size;i++){
        if(fixed.includes(i)) continue;
        fullU[i]=U[c++];
    }
    const reactions=multiplyMatrixVector(K,fullU).map((v,i)=>v-F[i]);
    return {nodes,displacements:fullU,reactions,supportIndices:supports};
}

function computeDiagrams(state,nodes,reactions){
    const events=[];
    const totalLength=nodes[nodes.length-1];
    for(let i=0;i<nodes.length;i++) events.push({x:nodes[i],P:reactions[2*i]});
    (state.pointLoads||[]).forEach(p=>events.push({x:p.x,P:p.P}));
    (state.lineLoads||[]).forEach(l=>events.push({x:l.start,w:l.w,start:true}));
    (state.lineLoads||[]).forEach(l=>events.push({x:l.end,w:l.w,end:true}));
    events.sort((a,b)=>a.x-b.x);

    let shear=0,moment=0;
    const xs=[0],shearVals=[0],momentVals=[0];
    let activeWs=[];
    for(let i=0;i<events.length;i++){
        const e=events[i];
        const dx=e.x-xs[xs.length-1];
        const wTotal=activeWs.reduce((a,b)=>a+b,0);
        const shearBefore=shear;
        shear-=wTotal*dx;
        moment+=shearBefore*dx -0.5*wTotal*dx*dx;
        xs.push(e.x); shearVals.push(shear); momentVals.push(moment);
        if(e.P){
            shear-=e.P;
            xs.push(e.x); shearVals.push(shear); momentVals.push(moment);
        }
        if(e.start) activeWs.push(e.w);
        if(e.end){
            const idx=activeWs.indexOf(e.w);
            if(idx>-1) activeWs.splice(idx,1);
        }
    }
    const lastX=xs[xs.length-1];
    if(lastX<totalLength){
        const dx=totalLength-lastX;
        const wTotal=activeWs.reduce((a,b)=>a+b,0);
        const shearBefore=shear;
        shear-=wTotal*dx;
        moment+=shearBefore*dx -0.5*wTotal*dx*dx;
        xs.push(totalLength); shearVals.push(shear); momentVals.push(moment);
    }
    return {xs,shearVals,momentVals};
}


if (typeof module !== 'undefined' && module.exports) {
    module.exports = { computeResults, computeDiagrams, setCrossSections, getSelfWeightLineLoads, getCrossSection, computeSectionDesign, computeInertia, computeWeakAxisInertia };
}

if (typeof window !== 'undefined') {
    window.computeResults = computeResults;
    window.computeDiagrams = computeDiagrams;
    window.setCrossSections = setCrossSections;
    window.getSelfWeightLineLoads = (spans,name)=>getSelfWeightLineLoads(spans,name);
    window.computeSectionDesign = computeSectionDesign;
    window.computeInertia = computeInertia;
    window.computeWeakAxisInertia = computeWeakAxisInertia;
}
