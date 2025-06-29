const assert = require('assert');
const {computeResults, computeSectionDesign} = require('../solver');

function close(actual, expected, tol, msg){
  if(Math.abs(actual-expected) > tol) throw new Error(msg+` expected ${expected} got ${actual}`);
}

(function testUniformLoad(){
  const L = 10;
  const w = 10; // N/m
  const E = 210e9;
  const I = 1e-6;
  const state = {spans:[{length:L}], pointLoads:[], lineLoads:[{w,start:0,end:L}], maxNodeDist:L/10, E, I};
  const {nodes,reactions,displacements} = computeResults(state);
  const Rexpected = -w*L/2;
  close(reactions[0], Rexpected, 1e-2, 'left reaction');
  close(reactions[2*10], Rexpected, 1e-2, 'right reaction');

  const midNode = nodes.findIndex(n=>Math.abs(n-L/2)<1e-6);
  const def = displacements[2*midNode];
  const defExpect = 5*w*Math.pow(L,4)/(384*E*I);
  close(def, defExpect, Math.abs(defExpect*0.05), 'midspan deflection');
})();

// Beam with left overhang (no left support)
(function testLeftOverhang(){
  const L = 10;
  const w = 10; // N/m
  const E = 210e9;
  const I = 1e-6;
  const state = {
    spans:[{length:L},{length:L}],
    pointLoads:[],
    lineLoads:[{w,start:0,end:2*L}],
    maxNodeDist:L/10,
    E,
    I,
    leftSupport:false,
    rightSupport:true
  };
  const {reactions} = computeResults(state);
  close(reactions[2*10], -w*2*L, 1e-2, 'interior support reaction');
  close(reactions[2*20], 0, 1e-6, 'end support reaction');
})();

// Steel design extra results
(function testSteelDesignExtras(){
  const design = computeSectionDesign('IPE100', {unbracedLength: 3});
  assert(design.Iw > 0, 'Iw not computed');
  assert(design.It > 0, 'It not computed');
  assert(design.Iz > 0, 'Iz not computed');
  assert(design.Mcr > 0, 'Mcr not computed');
  assert(design.chiLT > 0 && design.chiLT <= 1, 'chiLT invalid');
  assert(design.Lb === 3, 'Lb not stored');
  assert(design.E > 0 && design.G > 0, 'moduli not stored');
})();
