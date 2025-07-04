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

function invertMatrix(A){
    const n=A.length;
    const M=A.map(r=>r.slice());
    const I=Array.from({length:n},(_,i)=>{
        const row=new Array(n).fill(0); row[i]=1; return row;});
    for(let i=0;i<n;i++){
        let max=i;
        for(let k=i+1;k<n;k++) if(Math.abs(M[k][i])>Math.abs(M[max][i])) max=k;
        [M[i],M[max]]=[M[max],M[i]];
        [I[i],I[max]]=[I[max],I[i]];
        const pivot=M[i][i];
        for(let j=0;j<n;j++){M[i][j]/=pivot; I[i][j]/=pivot;}
        for(let k=0;k<n;k++) if(k!==i){
            const f=M[k][i];
            for(let j=0;j<n;j++){M[k][j]-=f*M[i][j]; I[k][j]-=f*I[i][j];}
        }
    }
    return I;
}

function frameElementWithReleases(E,A,I,L,rel){
    const kLocal=frameElementStiffness(E,A,I,L);
    const INF=1e12;
    const kx1=rel?.kx1===undefined?INF:(rel.kx1<0?INF:rel.kx1);
    const ky1=rel?.ky1===undefined?INF:(rel.ky1<0?INF:rel.ky1);
    const cz1=rel?.cz1===undefined?INF:(rel.cz1<0?INF:rel.cz1);
    const kx2=rel?.kx2===undefined?INF:(rel.kx2<0?INF:rel.kx2);
    const ky2=rel?.ky2===undefined?INF:(rel.ky2<0?INF:rel.ky2);
    const cz2=rel?.cz2===undefined?INF:(rel.cz2<0?INF:rel.cz2);
    const K=Array.from({length:12},()=>new Array(12).fill(0));
    for(let i=0;i<6;i++) for(let j=0;j<6;j++) K[6+i][6+j]=kLocal[i][j];
    const addSpr=(k,n,b)=>{
        if(k===0) return;
        K[n][n]+=k; K[b][b]+=k; K[n][b]-=k; K[b][n]-=k;
    };
    addSpr(kx1,0,6); addSpr(ky1,1,7); addSpr(cz1,2,8);
    addSpr(kx2,3,9); addSpr(ky2,4,10); addSpr(cz2,5,11);
    const Knn=K.slice(0,6).map(r=>r.slice(0,6));
    const Knb=K.slice(0,6).map(r=>r.slice(6));
    const Kbn=K.slice(6).map(r=>r.slice(0,6));
    const Kbb=K.slice(6).map(r=>r.slice(6));
    const KbbInv=invertMatrix(Kbb);
    const temp=multiplyMatrix(Knb,KbbInv);
    const sub=multiplyMatrix(temp,Kbn);
    const Kcond=Knn.map((row,i)=>row.map((v,j)=>v-sub[i][j]));
    return {Kcond,KbbInv,Kbn,kLocal};
}

function computeFrameResults(frame){
    const n=frame.nodes.length; if(n===0) return null;
    const dof=3*n;
    const K=Array.from({length:dof},()=>new Array(dof).fill(0));
    const F=new Array(dof).fill(0);
    frame.beams.forEach(el=>{
        if(el.on===false) return;
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy);
        const c=dx/L, s=dy/L;
        const E=el.E||frame.E||210e9;
        const I=el.I||frame.I||1e-6;
        const A=el.A||frame.A||0.001;
        const rel={kx1:el.kx1,ky1:el.ky1,cz1:el.cz1,kx2:el.kx2,ky2:el.ky2,cz2:el.cz2};
        const {Kcond}=frameElementWithReleases(E,A,I,L,rel);
        const T=[
            [ c, s,0, 0,0,0],
            [-s, c,0, 0,0,0],
            [ 0,0,1, 0,0,0],
            [ 0,0,0, c,s,0],
            [ 0,0,0,-s,c,0],
            [ 0,0,0, 0,0,1]
        ];
        const kG=multiplyMatrix(transpose(T),multiplyMatrix(Kcond,T));
        const dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2];
        for(let i=0;i<6;i++) for(let j=0;j<6;j++) K[dofs[i]][dofs[j]]+=kG[i][j];
    });
    frame.loads.forEach(l=>{
        if(l.Px) F[3*l.node]+=l.Px;
        if(l.Py) F[3*l.node+1]+=l.Py;
        if(l.Mz) F[3*l.node+2]+=l.Mz;
    });
    (frame.memberPointLoads||[]).forEach(l=>{
        const el=frame.beams[l.beam];
        if(!el||el.on===false) return;
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy); if(L===0) return;
        const c=dx/L, s=dy/L;
        const a=l.x; const b=L-a;
        const local=[0,0,0,0,0,0];
        const FxLocal=c*(l.Fx||0)+s*(l.Fy||0);
        const FyLocal=-s*(l.Fx||0)+c*(l.Fy||0);
        if(FxLocal){
            local[0]+=FxLocal*(1-a/L);
            local[3]+=FxLocal*(a/L);
        }
        if(FyLocal){
            const P=FyLocal;
            local[1]+=P*b*b*(3*a+b)/Math.pow(L,3);
            local[2]+=P*a*b*b/Math.pow(L,2);
            local[4]+=P*a*a*(3*b+a)/Math.pow(L,3);
            local[5]+=-P*a*a*b/Math.pow(L,2);
        }
        if(l.Mz){
            local[2]+=l.Mz*(1-a/L);
            local[5]+=l.Mz*(a/L);
        }
        const T=[[ c,-s,0,0,0,0],[ s, c,0,0,0,0],[0,0,1,0,0,0],[0,0,0, c,-s,0],[0,0,0, s, c,0],[0,0,0,0,0,1]];
        const gl=multiplyMatrixVector(transpose(T),local);
        const dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2];
        for(let i=0;i<6;i++) F[dofs[i]]+=gl[i];
    });
    (frame.memberLineLoads||[]).forEach(l=>{
        const el=frame.beams[l.beam];
        if(!el||el.on===false) return;
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy); if(L===0) return;
        const c=dx/L, s=dy/L;
        const local=[0,0,0,0,0,0];
        const isFullUniform = (l.start === 0 || l.start === undefined) &&
                              (l.end >= L || l.end === undefined) &&
                              l.wX1 === l.wX2 && l.wY1 === l.wY2;

        if (isFullUniform) {
            const wX = l.wX1 || 0;
            const wY = l.wY1 || 0;
            const w_ax = c * wX + s * wY;
            const w_prp = -s * wX + c * wY;

            if (w_ax !== 0) {
                local[0] += w_ax * L / 2;
                local[3] += w_ax * L / 2;
            }
            if (w_prp !== 0) {
                local[1] += w_prp * L / 2;
                local[2] += w_prp * L * L / 12;
                local[4] += w_prp * L / 2;
                local[5] += -w_prp * L * L / 12;
            }
        } else {
            const segs=50;
            const start=l.start || 0;
            const end=l.end || L;
            for(let i=0;i<segs;i++){
                const t1=i/segs, t2=(i+1)/segs;
                const x1=start + (end-start)*t1;
                const x2=start + (end-start)*t2;
                if(x2<=0||x1>=L) continue;
                const mid=(x1+x2)/2;
                const wX=l.wX1+(l.wX2-l.wX1)*((t1+t2)/2);
                const wY=l.wY1+(l.wY2-l.wY1)*((t1+t2)/2);
                const FxLocal=(c*wX + s*wY)*(x2-x1);
                const FyLocal=(-s*wX + c*wY)*(x2-x1);
                const a=mid; const b=L-a;
                if(FxLocal){
                    local[0]+=FxLocal*(1-a/L);
                    local[3]+=FxLocal*(a/L);
                }
                if(FyLocal){
                    const P=FyLocal;
                    local[1]+=P*b*b*(3*a+b)/Math.pow(L,3);
                    local[2]+=P*a*b*b/Math.pow(L,2);
                    local[4]+=P*a*a*(3*b+a)/Math.pow(L,3);
                    local[5]+=-P*a*a*b/Math.pow(L,2);
                }
            }
        }
        const T=[[ c,-s,0,0,0,0],[ s, c,0,0,0,0],[0,0,1,0,0,0],[0,0,0, c,-s,0],[0,0,0, s, c,0],[0,0,0,0,0,1]];
        const gl=multiplyMatrixVector(transpose(T),local);
        const dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2];
        for(let j=0;j<6;j++) F[dofs[j]]+=gl[j];
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
    const reactions=multiplyMatrixVector(K,full).map((v,i)=>v-F[i]);
    return {displacements:full,reactions};
}

function computeFrameDiagrams(frame,res,divisions=1){
    const diags=[];
    const segs=Math.max(1,10*divisions);
    frame.beams.forEach((el,idx)=>{
        if(el.on===false) return;
        const n1=el.n1,n2=el.n2;
        const p1=frame.nodes[n1], p2=frame.nodes[n2];
        const dx=p2.x-p1.x, dy=p2.y-p1.y;
        const L=Math.hypot(dx,dy);
        const c=dx/L, s=dy/L;
        const E=el.E||frame.E||210e9; const I=el.I||frame.I||1e-6; const A=el.A||frame.A||0.001;
        const rel={kx1:el.kx1,ky1:el.ky1,cz1:el.cz1,kx2:el.kx2,ky2:el.ky2,cz2:el.cz2};
        const {Kcond,KbbInv,Kbn,kLocal}=frameElementWithReleases(E,A,I,L,rel);
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
        const Ub=multiplyMatrixVector(KbbInv,multiplyMatrixVector(Kbn,dLocal)).map(v=>-v);
        const fLocal=multiplyMatrixVector(kLocal,Ub);

        const eq=[0,0,0,0,0,0];
        (frame.memberPointLoads||[]).filter(l=>l.beam===idx).forEach(l=>{
            const a=l.x; const b=L-a;
            const FxLocal=c*(l.Fx||0)+s*(l.Fy||0);
            const FyLocal=-s*(l.Fx||0)+c*(l.Fy||0);
            if(FxLocal){
                eq[0]+=FxLocal*(1-a/L);
                eq[3]+=FxLocal*(a/L);
            }
            if(FyLocal){
                const P=FyLocal;
                eq[1]+=P*b*b*(3*a+b)/Math.pow(L,3);
                eq[2]+=P*a*b*b/Math.pow(L,2);
                eq[4]+=P*a*a*(3*b+a)/Math.pow(L,3);
                eq[5]+=-P*a*a*b/Math.pow(L,2);
            }
            if(l.Mz){
                eq[2]+=l.Mz*(1-a/L);
                eq[5]+=l.Mz*(a/L);
            }
        });
        (frame.memberLineLoads||[]).filter(l=>l.beam===idx).forEach(l=>{
            const isFullUniform = (l.start === 0 || l.start === undefined) &&
                                  (l.end >= L || l.end === undefined) &&
                                  l.wX1 === l.wX2 && l.wY1 === l.wY2;

            if (isFullUniform) {
                const wX = l.wX1 || 0;
                const wY = l.wY1 || 0;
                const w_ax = c * wX + s * wY;
                const w_prp = -s * wX + c * wY;

                if (w_ax !== 0) {
                    eq[0] += w_ax * L / 2;
                    eq[3] += w_ax * L / 2;
                }
                if (w_prp !== 0) {
                    eq[1] += w_prp * L / 2;
                    eq[2] += w_prp * L * L / 12;
                    eq[4] += w_prp * L / 2;
                    eq[5] += -w_prp * L * L / 12;
                }
            } else {
                const segs=50;
                const start=l.start || 0;
                const end=l.end || L;
                for(let i=0;i<segs;i++){
                    const t1=i/segs, t2=(i+1)/segs;
                    const x1=start+(end-start)*t1;
                    const x2=start+(end-start)*t2;
                    if(x2<=0||x1>=L) continue;
                    const mid=(x1+x2)/2;
                    const wX=l.wX1+(l.wX2-l.wX1)*((t1+t2)/2);
                    const wY=l.wY1+(l.wY2-l.wY1)*((t1+t2)/2);
                    const FxLocal=(c*wX + s*wY)*(x2-x1);
                    const FyLocal=(-s*wX + c*wY)*(x2-x1);
                    const a=mid; const b=L-a;
                    if(FxLocal){
                        eq[0]+=FxLocal*(1-a/L);
                        eq[3]+=FxLocal*(a/L);
                    }
                    if(FyLocal){
                        const P=FyLocal;
                        eq[1]+=P*b*b*(3*a+b)/Math.pow(L,3);
                        eq[2]+=P*a*b*b/Math.pow(L,2);
                        eq[4]+=P*a*a*(3*b+a)/Math.pow(L,3);
                        eq[5]+=-P*a*a*b/Math.pow(L,2);
                    }
                }
            }
        });

        const fAdj=fLocal.map((v,i)=>v+eq[i]);
        if(rel.kx1===0) fAdj[0]=0;
        if(rel.ky1===0) fAdj[1]=0;
        if(rel.cz1===0) fAdj[2]=0;
        if(rel.kx2===0) fAdj[3]=0;
        if(rel.ky2===0) fAdj[4]=0;
        if(rel.cz2===0) fAdj[5]=0;
        // Use a standard sign convention where positive internal forces
        // correspond to positive diagrams without extra sign inversions.
        // Right-end forces retain the same sign convention as the left end.
        const N1=fAdj[0], N2=fAdj[3];
        const V1=fAdj[1], V2=fAdj[4];
        const M1=fAdj[2], M2=fAdj[5];

        const pointLoads=[{x:L,Fx:-N2,Fy:-V2,M:-M2}];
        (frame.memberPointLoads||[]).filter(l=>l.beam===idx).forEach(l=>{
            const FxLocal=c*(l.Fx||0)+s*(l.Fy||0);
            const FyLocal=-s*(l.Fx||0)+c*(l.Fy||0);
            pointLoads.push({x:l.x,Fx:FxLocal,Fy:FyLocal,M:l.Mz||0});
        });
        const lineSegs=[];
        (frame.memberLineLoads||[]).filter(l=>l.beam===idx).forEach(l=>{
            for(let i=0;i<segs;i++){
                const t1=i/segs, t2=(i+1)/segs;
                const x1=l.start+(l.end-l.start)*t1;
                const x2=l.start+(l.end-l.start)*t2;
                if(x2<=0||x1>=L) continue;
                const wX=l.wX1+(l.wX2-l.wX1)*((t1+t2)/2);
                const wY=l.wY1+(l.wY2-l.wY1)*((t1+t2)/2);
                const FxLocal=c*wX + s*wY;
                const FyLocal=-s*wX + c*wY;
                lineSegs.push({start:x1,end:x2,wX:FxLocal,wY:FyLocal});
            }
        });

        const positions=new Set([0,L]);
        for(let d=1; d<divisions; d++) positions.add(L*d/divisions);
        pointLoads.forEach(p=>positions.add(p.x));
        lineSegs.forEach(seg=>{positions.add(seg.start); positions.add(seg.end);});
        const posArr=Array.from(positions).sort((a,b)=>a-b);
        // Start diagrams with the actual member end forces
        let shear=V1, moment=M1, normal=N1;
        const shearArr=[{x:0,y:shear}], momentArr=[{x:0,y:moment}], normalArr=[{x:0,y:normal}];
        let active=[];
        for(let i=0;i<posArr.length-1;i++){
            const x1=posArr[i]; const x2=posArr[i+1];
            const dx=x2-x1;
            const wX=active.reduce((a,b)=>a+b.wX,0);
            const wY=active.reduce((a,b)=>a+b.wY,0);
            const shearBefore=shear; const normalBefore=normal;
            shear-=wY*dx;
            normal-=wX*dx;
            moment+=shearBefore*dx-0.5*wY*dx*dx;
            shearArr.push({x:x2,y:shear});
            momentArr.push({x:x2,y:moment});
            normalArr.push({x:x2,y:normal});
            pointLoads.filter(p=>Math.abs(p.x-x2)<1e-8).forEach(p=>{
                shear-=p.Fy; normal-=p.Fx; moment-=p.M; shearArr.push({x:x2,y:shear}); momentArr.push({x:x2,y:moment}); normalArr.push({x:x2,y:normal});
            });
            lineSegs.filter(seg=>Math.abs(seg.start-x2)<1e-8).forEach(seg=>active.push(seg));
            lineSegs.filter(seg=>Math.abs(seg.end-x2)<1e-8).forEach(seg=>{const idx=active.indexOf(seg); if(idx>-1) active.splice(idx,1);});
        }
        diags.push({shear:shearArr,moment:momentArr,normal:normalArr});
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
