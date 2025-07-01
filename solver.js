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

function transpose(M){
    const r=[]; for(let i=0;i<M[0].length;i++){r[i]=[]; for(let j=0;j<M.length;j++) r[i][j]=M[j][i];}
    return r;
}

function multiplyMatrix(A,B){
    const r=[]; for(let i=0;i<A.length;i++){r[i]=[]; for(let j=0;j<B[0].length;j++){let sum=0; for(let k=0;k<B.length;k++) sum+=A[i][k]*B[k][j]; r[i][j]=sum;}}
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
        if(cs.I_y !== undefined) return cs.I_y;
        const b = cs.b_mm/1000;
        const h = cs.h_mm/1000;
        return b*Math.pow(h,3)/12;
    }
    const h = cs.h_mm !== undefined ? cs.h_mm/1000 : undefined;
    const b = cs.b_mm !== undefined ? cs.b_mm/1000 : undefined;
    const tw = cs.tw_mm !== undefined ? cs.tw_mm/1000 : undefined;
    const tf = cs.tf_mm !== undefined ? cs.tf_mm/1000 : undefined;
    if(h===undefined || b===undefined || tw===undefined || tf===undefined)
        return cs.I_y !== undefined ? cs.I_y : (cs.Iy_m4 !== undefined ? cs.Iy_m4 : (cs.Iz_m4 !== undefined ? cs.Iz_m4 : (cs.I_z !== undefined ? cs.I_z : 0)));
    const hw = h - 2*tf;
    const Iweb = tw*Math.pow(hw,3)/12;
    const If = b*Math.pow(tf,3)/12;
    const d = hw/2 + tf/2;
    const Icalc = 2*(If + b*tf*d*d) + Iweb; // dimensions already in m
    return Icalc;
}

function computeWeakAxisInertia(cs){
    if(!cs) return 0;
    if(cs.Iz_m4 !== undefined) return cs.Iz_m4;
    if(cs.I_z !== undefined) return cs.I_z;
    if(cs.series === 'sawn' || cs.series === 'glulam'){
        const b = cs.b_mm/1000;
        const h = cs.h_mm/1000;
        return h*Math.pow(b,3)/12;
    }
    const h = cs.h_mm !== undefined ? cs.h_mm/1000 : undefined;
    const b = cs.b_mm !== undefined ? cs.b_mm/1000 : undefined;
    const tw = cs.tw_mm !== undefined ? cs.tw_mm/1000 : undefined;
    const tf = cs.tf_mm !== undefined ? cs.tf_mm/1000 : undefined;
    if(h===undefined || b===undefined || tw===undefined || tf===undefined)
        return 0;
    const hw = h - 2*tf;
    const Iweb = hw*Math.pow(tw,3)/12;
    const If = b*Math.pow(tf,3)/12;
    return 2*If + Iweb;
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
    const Iz = computeWeakAxisInertia(cs);
    let Iw = (tf * Math.pow(b, 3) / 24) * Math.pow(h - tf, 2);
    let Mcr, chiLT;
    let Lb, G = opts.G !== undefined ? opts.G : 81e9;
    let C1 = opts.C1 !== undefined ? opts.C1 : 1.0;
    let C2 = opts.C2 !== undefined ? opts.C2 : 0.0;
    let C3 = opts.C3 !== undefined ? opts.C3 : 0.0;
    let kw = opts.kw !== undefined ? opts.kw : 1.0;
    if(typeof opts.unbracedLength === 'number' && opts.unbracedLength > 0){
        Lb = opts.unbracedLength;
        const term_sqrt_1 = (Math.PI * Math.PI * E * Iz) / (Lb * Lb);
        const term_sqrt_2 = Iw / Iz;
        const term_sqrt_3 = (Lb * Lb * G * It) / (Math.PI * Math.PI * E * Iz);
        Mcr = C1 * term_sqrt_1 * Math.sqrt(term_sqrt_2 + term_sqrt_3);
        const lambdaRel = Math.sqrt((fy*W/gammaM0)/Mcr);
        const alpha = 0.34;
        const phi = 0.5*(1 + alpha*(lambdaRel-0.2) + lambdaRel*lambdaRel);
        chiLT = 1/(phi + Math.sqrt(phi*phi - lambdaRel*lambdaRel));
        MRdLBA = chiLT*(fy*W/gammaM0);
    }
    return {EI, MRd, MRdLBA, VRd, W, gamma: gammaM0, material: 'steel',
            Iw, It, Iz, Lb, E, G, C1, C2, C3, kw, Mcr, chiLT};
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



function frameElementStiffness(E,A,I,L){
    const a=E*A/L;
    const b=12*E*I/Math.pow(L,3);
    const c=6*E*I/Math.pow(L,2);
    const d=4*E*I/L;
    const e=2*E*I/L;
    return [
        [ a, 0, 0, -a, 0, 0],
        [ 0, b, c, 0,-b, c],
        [ 0, c, d, 0,-c, e],
        [-a, 0, 0,  a, 0, 0],
        [ 0,-b,-c,0, b,-c],
        [ 0, c, e, 0,-c, d]
    ];
}

function computeFrameResults(frame){
    const n=frame.nodes.length; if(n===0) return null;
    const dof=3*n;
    const K=Array.from({length:dof},()=>new Array(dof).fill(0));
    const F=new Array(dof).fill(0);
    frame.beams.forEach(el=>{
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy);
        const c=dx/L, s=dy/L;
        const E=el.E||frame.E||210e9;
        const I=el.I||frame.I||1e-6;
        const A=el.A||frame.A||0.001;
        const kLocal=frameElementStiffness(E,A,I,L);
        const T=[
            [ c, s,0, 0,0,0],
            [-s, c,0, 0,0,0],
            [ 0,0,1, 0,0,0],
            [ 0,0,0, c,s,0],
            [ 0,0,0,-s,c,0],
            [ 0,0,0, 0,0,1]
        ];
        const kG=multiplyMatrix(transpose(T),multiplyMatrix(kLocal,T));
        const dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2];
        for(let i=0;i<6;i++) for(let j=0;j<6;j++) K[dofs[i]][dofs[j]]+=kG[i][j];
    });
    frame.loads.forEach(l=>{
        if(l.Px) F[3*l.node]+=l.Px;
        if(l.Py) F[3*l.node+1]+=l.Py;
        if(l.Mz) F[3*l.node+2]+=l.Mz;
    });
    const fixed=[];
    frame.supports.forEach(s=>{
        if(s.fixX) fixed.push(3*s.node);
        if(s.fixY) fixed.push(3*s.node+1);
        if(s.fixRot) fixed.push(3*s.node+2);
    });
    const {Kmod,Fmod}=applyBC(K,F,fixed);
    const U=gaussSolve(Kmod,Fmod);
    const full=new Array(dof).fill(0); let c=0;
    for(let i=0;i<dof;i++){if(!fixed.includes(i)) full[i]=U[c++];}
    return {displacements:full};
}

function computeFrameDiagrams(frame,res,divisions=1){
    const diags=[];
    frame.beams.forEach(el=>{
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy);
        const c=dx/L, s=dy/L;
        const E=el.E||frame.E||210e9; const I=el.I||frame.I||1e-6; const A=el.A||frame.A||0.001;
        const kLocal=frameElementStiffness(E,A,I,L);
        const T=[
            [ c, s,0, 0,0,0],
            [-s, c,0, 0,0,0],
            [ 0,0,1, 0,0,0],
            [ 0,0,0, c,s,0],
            [ 0,0,0,-s,c,0],
            [ 0,0,0, 0,0,1]
        ];
        const dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2];
        const dGlobal=dofs.map(i=>res.displacements[i]);
        const dLocal=multiplyMatrixVector(T,dGlobal);
        const fLocal=multiplyMatrixVector(kLocal,dLocal);
        const N1=-fLocal[0], N2=-fLocal[3];
        const V1=fLocal[1], V2=-fLocal[4];
        const M1=fLocal[2], M2=-fLocal[5];
        const shear=[], moment=[], normal=[];
        for(let j=0;j<=divisions;j++){
            const t=j/divisions;
            shear.push({x:t*L,y:V1+(V2-V1)*t});
            moment.push({x:t*L,y:M1+(M2-M1)*t});
            normal.push({x:t*L,y:N1+(N2-N1)*t});
        }
        diags.push({shear,moment,normal});
    });
    return diags;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { computeResults, computeDiagrams, computeFrameResults, computeFrameDiagrams, setCrossSections, getSelfWeightLineLoads, getCrossSection, computeSectionDesign, computeInertia, computeWeakAxisInertia };
}

if (typeof window !== 'undefined') {
    window.computeResults = computeResults;
    window.computeDiagrams = computeDiagrams;
    window.computeFrameResults = computeFrameResults;
    window.computeFrameDiagrams = computeFrameDiagrams;
    window.setCrossSections = setCrossSections;
    window.getSelfWeightLineLoads = (spans,name)=>getSelfWeightLineLoads(spans,name);
    window.computeSectionDesign = computeSectionDesign;
    window.computeInertia = computeInertia;
    window.computeWeakAxisInertia = computeWeakAxisInertia;
}
