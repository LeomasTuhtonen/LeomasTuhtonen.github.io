const assert = require('assert');
const {computeResults, computeSectionDesign, computeFrameResults, computeFrameResultsPDelta_Kg, computeFrameResultsLBA, computeFrameDiagrams} = require('../solver');

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

// Concrete design basic test
(function testConcreteDesign(){
  const design = computeSectionDesign('200x400', {material:'concrete', rho:0.01});
  assert(design.MRd > 0, 'MRd not computed');
  assert(design.VRd > 0, 'VRd not computed');
})();

// Simple frame test
(function testSimpleFrame(){
  const frame = {
    nodes:[{x:0,y:0},{x:0,y:3},{x:4,y:3},{x:4,y:0}],
    beams:[{n1:0,n2:1},{n1:1,n2:2},{n1:2,n2:3}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true},{node:3,fixX:true,fixY:true,fixRot:true}],
    loads:[{node:1,Py:-1000}]
  };
  const res = computeFrameResults(frame);
  assert(res.displacements.length === frame.nodes.length*3, 'frame result size');
})();

// Frame member load test
(function testFrameMemberPoint(){
  const frame={
    nodes:[{x:0,y:0},{x:0,y:3}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true}],
    loads:[],
    memberPointLoads:[{beam:0,x:1,Fy:-1000}]
  };
  const res=computeFrameResults(frame);
  assert(res.displacements.length===frame.nodes.length*3,'member load result');
})();

// Global axis point load transformation
(function testGlobalAxisPointLoad(){
  const frame={
    nodes:[{x:0,y:0},{x:1,y:0}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true},{node:1,fixY:true,fixRot:true}],
    loads:[],
    memberPointLoads:[{beam:0,x:0.5,Fx:1000}]
  };
  const res=computeFrameResults(frame);
  assert(res.displacements[3]>0,'positive axial displacement expected');
})();

// Frame diagrams with member point load
(function testFrameDiagramsPoint(){
  const frame={
    nodes:[{x:0,y:0},{x:1,y:0}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true}],
    loads:[],
    memberPointLoads:[{beam:0,x:0.5,Fy:-1}]
  };
  const res=computeFrameResults(frame);
  const diags=computeFrameDiagrams(frame,res,1);
  const shear=diags[0].shear.map(p=>p.y);
  const moment=diags[0].moment.map(p=>p.y);
  const lastShear=shear[shear.length-1];
  const lastMoment=moment[moment.length-1];
  assert(Math.abs(shear[0]+1)<1e-4 && Math.abs(lastShear-0)<1e-4,'shear diagram');
  assert(Math.abs(moment[0]-0.5)<1e-4 && Math.abs(lastMoment-0)<1e-4,'moment diagram');
})();

// Moment release at beam start
(function testMomentRelease(){
  const frame={
    nodes:[{x:0,y:0},{x:1,y:0}],
    beams:[{n1:0,n2:1,cz1:0}],
    supports:[{node:0,fixX:true,fixY:true},{node:1,fixY:true}],
    loads:[],
    memberPointLoads:[{beam:0,x:0.5,Fy:-1}]
  };
  const res=computeFrameResults(frame);
  const diags=computeFrameDiagrams(frame,res,1);
  const startMoment=diags[0].moment[0].y;
  assert(Math.abs(startMoment-0) < 1e-6, 'moment value mismatch');
})();

// Uniform load on inclined member
(function testInclinedUniform(){
  const L = Math.sqrt(2);
  const frame={
    nodes:[{x:0,y:0},{x:1,y:1}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true},{node:1,fixX:true,fixY:true,fixRot:true}],
    loads:[],
    memberLineLoads:[{beam:0,start:0,end:L,wY1:-1,wY2:-1}]
  };
  const res=computeFrameResults(frame);
  close(res.reactions[1], L/2, 1e-6, 'left vertical reaction');
  close(res.reactions[4], L/2, 1e-6, 'right vertical reaction');
})();

// P-Delta effect increases displacement
(function testPDeltaDisplacement(){
  const frame={
    nodes:[{x:0,y:0},{x:0,y:3}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true}],
    loads:[{node:1,Px:1000,Py:-10000}]
  };
  const first=computeFrameResults(frame);
  const second=computeFrameResultsPDelta_Kg(frame);
  assert(Math.abs(second.displacements[3])>Math.abs(first.displacements[3]),'P-Delta should increase sway');
  const diag=second.diagrams || computeFrameDiagrams(frame,second,1);
  const nonZero=diag[0].moment.some(pt=>Math.abs(pt.y)>1e-8);
  assert(nonZero,'P-Delta diagrams should not be zero');
})();

// Basic LBA test
(function testLBA(){
  const frame={
    nodes:[{x:0,y:0},{x:0,y:1}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true},{node:1,fixX:true}],
    loads:[{node:1,Py:-1000}]
  };
  const res=computeFrameResultsLBA(frame,0);
  assert(res.alpha>1000,'alpha magnitude');
  assert(res.displacements.length===6,'mode displacement size');
  for(let i=1;i<res.modes.length;i++){
    assert(res.modes[i]>=res.modes[i-1],'modes sorted');
  }
})();

// Point load reaction should equal applied load
(function testCantileverPointLoad(){
  const frame={
    nodes:[{x:0,y:0},{x:1,y:0}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixX:true,fixY:true,fixRot:true}],
    loads:[],
    memberPointLoads:[{beam:0,x:0.5,Fy:-1000}]
  };
  const res=computeFrameResults(frame);
  close(res.reactions[1], 1000, 1e-6, 'cantilever vertical reaction');
})();

// Uniform line load symmetry
(function testUniformLineLoadSymmetry(){
  const L=2;
  const frame={
    nodes:[{x:0,y:0},{x:L,y:0}],
    beams:[{n1:0,n2:1}],
    supports:[{node:0,fixY:true,fixRot:true},{node:1,fixY:true,fixRot:true}],
    loads:[],
    memberLineLoads:[{beam:0,start:0,end:L,wY1:-1000,wY2:-1000}]
  };
  const res=computeFrameResults(frame);
  close(res.reactions[1], 1000*L/2, 1e-6, 'left uniform reaction');
  close(res.reactions[4], 1000*L/2, 1e-6, 'right uniform reaction');
})();
