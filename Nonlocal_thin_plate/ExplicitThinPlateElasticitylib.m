(* ::Package:: *)

Needs["NDSolve`FEM`"]
iDf={1, 10, 2, 11, 20};
iDfsum={1, 1, 2, 2, 2};
iDffact={1, 1, 2, 1, 2};
orderDeri=2;
numDeri=5;
Py={y, x, y^2, x*y, x^2};
AbortOnMessage[t_:True] := If[t, messageHandler = If[Last[#1], Abort[]] & ; 
     Internal`AddHandler["Message", messageHandler], 
    Internal`RemoveHandler["Message", messageHandler]]
BCValuesFunc[pL_List, coord_List, fun_, udim_] := 
   Module[{v, ndim, fdim}, ndim = Length[coord[[1]]]; 
     fdim = fun[coord[[pL[[1]]]]]; If[udim != fdim, 
      Print["Error, the dimension of fun must be the same as the number of \
multiple fields"]; Return[]]; v = ConstantArray[0., {Length[pL], udim}]; 
     Do[v[[i]] += fun[coord[[pL[[i]]]]], {i, Length[pL]}]; v]
DofIndexValue[pL_List, vl_List] := Module[{udim, pi, pu, dofIndex, dofValue, 
     k = 1}, udim = Length[vl[[1]]]; dofIndex = ConstantArray[0, 
       udim*Length[pL]]; dofValue = ConstantArray[0., udim*Length[pL]]; 
     pi = Range[udim]; Do[pu = pi + udim*(pL[[i]] - 1); 
       Do[If[vl[[i,j]] < 10^10, dofIndex[[k]] = pu[[j]]; dofValue[[k++]] = 
           vl[[i,j]]], {j, udim}], {i, Length[pL]}]; 
     {dofIndex[[1 ;; k - 1]], dofValue[[1 ;; k - 1]]}]
FindPoints[v1_List, xmin_List, xmax_List, range_] := 
   Module[{ilist = {}, y, r, ndim = Length[v1[[1]]], len = Length[v1]}, 
    Do[If[LessThan2[v1[[i]], xmax] && LessThan2[xmin, v1[[i]]], 
       AppendTo[ilist, i]], {i, range}]; ilist]
 
FindPoints[v1_List, xmin_List, xmax_List, show_:False] := 
   Module[{ilist = {}, y, r, ndim = Length[v1[[1]]], len = Length[v1]}, 
    Do[If[LessThan2[v1[[i]], xmax] && LessThan2[xmin, v1[[i]]], 
       AppendTo[ilist, i]], {i, len}]; 
     If[show, If[Length[v1[[1]]] == 2, Print[ListPlot[v1[[ilist]], 
          DataRange -> Automatic]]]; If[Length[v1[[1]]] == 3, 
        Print[ListPointPlot3D[v1[[ilist]]], DataRange -> Full]]]; ilist]
 
FindPoints[v1_List, fun_, dist_, show_:False] := 
   Module[{ilist = {}, y, r, ndim = Length[v1[[1]]], len = Length[v1]}, 
    Do[If[Abs[fun[v1[[i]]]] <= dist, AppendTo[ilist, i]], {i, len}]; 
     If[show, If[Length[v1[[1]]] == 2, Print[ListPlot[v1[[ilist]]]]]; 
       If[Length[v1[[1]]] == 3, Print[ListPointPlot3D[v1[[ilist]]]]]]; ilist]
LessThan2[v1_List, v2_List] := Module[{p = True}, 
    Do[If[v1[[i]] > v2[[i]], p = False; Break[]], {i, Length[v1]}]; Return[p]]
fixFun[xyz_List] := ConstantArray[0., udim]
forceFun[xyz_List] := {fx*dx^2}
GridNdim[num_, ndim_] := 
   Module[{idx = {"i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", 
       "t"}, id2 = {}, st = "Table[", al}, idx = idx[[1 ;; ndim]]; 
     st = StringJoin[st, ToString[idx]]; 
     Do[id2 = StringJoin[id2, ",", ToString[{idx[[i]], 0, num}]], {i, ndim}]; 
     al = Flatten[ToExpression[StringJoin[st, id2, "]"]]]; 
     ArrayReshape[al, {Length[al]/ndim, ndim}]]
 
GridNdim[xmin_List, xmax_List, dx_] := 
   Module[{idx = {"i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", 
       "t"}, id2 = {}, st = "Table[", al, nx = Ceiling[(xmax - xmin)/dx], 
     ndx, dxL, ndim, res}, ndim = Length[xmin]; idx = idx[[1 ;; ndim]]; 
     ndx = (xmax - xmin)/nx; st = StringJoin[st, ToString[idx]]; 
     Do[id2 = StringJoin[id2, ",", ToString[{idx[[i]], 0, nx[[i]]}]], 
      {i, ndim}]; al = Flatten[ToExpression[StringJoin[st, id2, "]"]]]; 
     res = ArrayReshape[al, {Length[al]/ndim, ndim}]; 
     Do[res[[All,i]] *= ndx[[i]], {i, ndim}]; Do[res[[i]] += xmin, 
      {i, Length[res]}]; res]
MyDelaunayMesh[coord_] := Block[{mm, mm2, ndim, ige}, 
    ndim = Length[coord[[1]]]; If[ndim < 2 || ndim > 3, 
      Print["Error, delaunayMesh only applies for 2D and 3D"]; Abort[]; ]; 
     mm = DelaunayMesh[coord]; mm2 = If[ndim == 2, Cases[MeshCells[mm, 2], 
        Polygon[x_] :> x, -1], Cases[MeshCells[mm, 3], Tetrahedron[x_] :> x, 
        -1]]; ige = SimplexMeshCheckQuality[coord, mm2]; 
     If[ige[[1]] != -1, mm2 = mm2[[ige]]]; Abaqus2Mesh[coord, mm2]]
 
MyDelaunayMesh[coord_, uvw_] := Block[{mm, mm2, ndim}, 
    mm = DelaunayMesh[coord]; mm2 = If[ndim == 2, Cases[MeshCells[mm, 2], 
        Polygon[x_] :> x, -1], Cases[MeshCells[mm, 3], Tetrahedron[x_] :> x, 
        -1]]; Abaqus2Mesh[coord + uvw, mm2]]
SimplexMeshCheckQuality[coord_, ele_] := 
   Block[{ndim = Length[coord[[1]]], ne = Length[ele], xy, vl, ei, vm = 0., 
     ige}, vl = Table[0., {ne}]; If[ndim == 2, 
      Do[ei = ele[[i]]; xy = coord[[ei]]; vl[[i]] = AreaOfTriangle[xy]; , 
       {i, ne}]]; If[ndim == 3, Do[ei = ele[[i]]; xy = coord[[ei]]; 
        vl[[i]] = VolumeOfTet[xy]; , {i, ne}]]; vm = Mean[vl]/10^3; 
     ige = Flatten[Position[vl, _?(#1 > vm & )]]; If[Length[ige] < ne, 
      Print[ne - Length[ige], " bad elements were detected and deleted!"], 
      Print["All elements are in good condition"]; ige = {-1}; ]; ige]
AreaOfTriangle[v_] := Block[{v12 = v[[2]] - v[[1]], v13 = v[[3]] - v[[1]]}, 
    0.5*Abs[v12[[1]]*v13[[2]] - v12[[2]]*v13[[1]]]]
VolumeOfTet[v_] := Block[{a11, a12, a13, a21, a22, a23, a31, a32, a33}, 
    {a11, a12, a13} = v[[4]] - v[[1]]; {a21, a22, a23} = v[[3]] - v[[1]]; 
     {a31, a32, a33} = v[[2]] - v[[1]]; 
     (1./6)*Abs[(-a13)*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - 
        a12*a21*a33 + a11*a22*a33]]
Abaqus2Mesh[coord_, mesh_] := Block[{ec, el = {}}, 
    ec = ElementClassify[mesh, Length[coord[[1]]]]; 
     If[Length[ec[[1]]] > 0, AppendTo[el, TriangleElement[mesh[[ec[[1]]]]]]]; 
     If[Length[ec[[2]]] > 0, AppendTo[el, QuadElement[mesh[[ec[[2]]]]]]]; 
     If[Length[ec[[3]]] > 0, AppendTo[el, TetrahedronElement[
        mesh[[ec[[3]]]]]]]; If[Length[ec[[4]]] > 0, 
      AppendTo[el, HexahedronElement[mesh[[ec[[4]]]]]]]; 
     ToElementMesh["Coordinates" -> coord, "MeshElements" -> el]]
ElementClassify[mesh_, ndim_] := Block[{pos = ConstantArray[0, Length[mesh]], 
     toqh = {{}, {}, {}, {}}}, 
    Do[pos[[i]] = ElementClassifyEc[Length[mesh[[i]]], ndim]; , 
      {i, Length[mesh]}]; toqh[[1]] = Flatten[Position[pos, 1]]; 
     toqh[[2]] = Flatten[Position[pos, 2]]; 
     toqh[[3]] = Flatten[Position[pos, 3]]; 
     toqh[[4]] = Flatten[Position[pos, 4]]; toqh]
ElementClassifyEc[len_, n_] := If[n == 2, If[len == 3 || len == 6, 1, 2], 
    If[len == 4 || len == 10, 3, 4]]
nearestNeighbours[coord_List, numNei_] := Block[{cfun, nei}, 
    If[$VersionNumber > 10, Nearest[coord -> Automatic, coord, numNei + 1][[
      All,2 ;; All]], cfun = Nearest[coord -> Automatic]; 
      Table[cfun[coord[[i]], numNei + 1], {i, Length[coord]}][[All,
       2 ;; All]]]]
 
nearestNeighbours[coord_List, nodes_List, numNei_] := 
   Block[{c1 = coord[[nodes]], n1}, n1 = nearestNeighbours[c1, numNei]; 
     n1 //. ListRule[Range[nodes], nodes]]
 
nearestNeighbours[coord_List, numNei_, dx_List] := 
   Block[{pl, dxMax = Max[dx], pi, pi2, dpi}, 
    pl = Nearest[coord -> Automatic, coord, {Infinity, 1.5*dxMax}]; 
     Do[If[Round[dxMax/dx[[i]]] == 1, Continue[]]; pi = {}; pi2 = pl[[i]]; 
       Do[dpi = Norm[coord[[pi2[[j]]]] - coord[[i]]]; If[dpi < 1.5*dx[[i]], 
          AppendTo[pi, pi2[[j]]]]; , {j, 2, Length[pi2]}]; pl[[i]] = pi; , 
      {i, Length[coord]}]; pl]
ListRule[var_List, values_List] := Apply[Rule, Transpose[{var, values}], {1}]
NOMPartialU0[uvw0_List, PW_List, Nei_List] := 
   Block[{isIinclude, iuvw, udim = Length[uvw0[[1]]], len = Length[uvw0], 
     res, uvw}, isIinclude = If[Nei[[1,1]] == 1, True, False]; 
     If[Length[Dimensions[uvw0]] == 1, uvw = Table[{i}, {i, uvw0}], 
      uvw = uvw0]; If[isIinclude, iuvw = Table[uvw[[Nei[[i]]]], {i, len}], 
      iuvw = Table[uvw[[Prepend[Nei[[i]], i]]], {i, len}]]; 
     Print[Dimensions[PW], " ", Dimensions[iuvw]]; 
     res = ListMatrixMultiply[PW, iuvw]; If[udim > 1, 
      res = Table[Transpose[res[[i]]], {i, len}]]; res]
NOMPwKhg[np0_, coord_List, vol_, Nei_List, pfun_, WeiF_, hgPen_] := 
   Module[{tl, ndof, Nnode, ndim, udim, kk, np, tl2}, 
    If[np0 < 0, np = Ceiling[Length[coord]/1000], np = np0]; 
     tl2 = RandomInteger[Length[coord], 10]; 
     tl2 = NOMPwKhgPartConditionNumber[tl2, coord, vol, Nei, pfun, WeiF, 
       hgPen]; Print["Condition number of NOM:", tl2]; 
     If[Max[tl2] > 10^12, Print["Warning, condition number too high, please \
consider more particles or weigthing function"]]; 
     tl = Ngroup[Length[coord], np]; 
     kk = ParallelTable[NOMPwKhgPart[tl[[i]], coord, vol, Nei, pfun, WeiF, 
        hgPen], {i, np}]; {Flatten[kk[[All,1]], 1], Flatten[kk[[All,2]], 1]}]
QuinticKer[r0_, h_] := With[{r = r0/h}, If[1 >= r, (1 - r)^5, 0] - 
     6*If[2./3 >= r, (2./3 - r)^5, 0] + 15*If[1./3 >= r, (1./3 - r)^5, 0]]
QuinticKer[r0_, h_] := With[{r = r0/h}, If[1 >= r, (1 - r)^5, 0] - 
     6*If[2./3 >= r, (2./3 - r)^5, 0] + 15*If[1./3 >= r, (1./3 - r)^5, 0]]
NOMPwKhgPartConditionNumber[nodeList_, coord_List, vol_, Nei_List, pfun_, 
    WeiF_, hgPen_] := Module[{pvol, Nnode, pwList, hgList, kp = 1, coordi, 
     pw, khg, NeiI, hg, hi, voli = 0, cL = {}}, Nnode = Length[coord]; 
     If[Length[vol] == 0, pvol = ConstantArray[vol, Nnode]; ]; 
     If[Length[vol] == Nnode, pvol = vol]; 
     pwList = EmptyList[Length[nodeList]]; 
     hgList = EmptyList[Length[nodeList]]; 
     Do[NeiI = Nei[[i]]; coordi = ListMinus[coord[[NeiI]], coord[[i]]]; 
       voli = pvol[[NeiI]]; hi = Norm[coordi[[-1]]]; 
       AppendTo[cL, KhgsConditionNumber[coordi, pfun, WeiF, voli, hgPen, 
         hi]]; , {i, nodeList}]; Print["For selected nodes:", nodeList]; cL]
EmptyList[n_Integer] := Block[{el = {}}, Do[AppendTo[el, {}], {i, n}]; el]
 
EmptyList[n_Integer, m_Integer] := Block[{el = {}}, 
    Do[AppendTo[el, {}]; Do[AppendTo[el[[i]], {}], {j, m}], {i, n}]; el]
ListMinus[v1_List, v0_List] := Module[{v2 = v1}, 
    Do[v2[[i]] -= v0, {i, Length[v1]}]; v2]
KhgsConditionNumber[v1_List, pfun_, wfun_, vol_List, penalty_, h_] := 
   Module[{len = Length[pfun[v1[[1]], 1.]], num = Length[v1], k, p, trH = 0., 
     r1, p1, w1, pkp, sc, ssc, khg, pw0, pw, wl, cll = {}}, 
    pkp = ConstantArray[0., {num, num}]; k = ConstantArray[0., {len, len}]; 
     p = ConstantArray[0., {num, len}]; 
     khg = ConstantArray[0., {num + 1, num + 1}]; wl = Norm /@ v1; 
     Do[wl[[i]] = wfun[wl[[i]], 1.1*h]*vol[[i]], {i, num}]; 
     wl *= 1./Total[wl]; Do[r1 = Norm[v1[[i]]]; w1 = wl[[i]]; 
       trH += w1*r1*r1; pkp[[i,i]] = w1; p1 = pfun[v1[[i]], h]; 
       p[[i]] = w1*p1; k += w1*p1 \[TensorProduct] p1, {i, num}]; 
     Eigenvalues[k, 1][[1]]/Eigenvalues[k, -1][[1]]]
Ngroup[nmax_, ngroup_] := Module[{tl = {}, i1, i2, inc}, 
    inc = Round[nmax/ngroup]; Do[AppendTo[tl, Range[inc*(i - 1) + 1, inc*i]], 
      {i, 1, ngroup - 1}]; AppendTo[tl, Range[inc*(ngroup - 1) + 1, nmax]]; 
     tl]
NOMPwKhgPart[nodeList_, coord_List, vol_, Nei_List, pfun_, WeiF_, hgPen_] := 
   Module[{pvol, Nnode, pwList, hgList, kp = 1, coordi, pw, khg, NeiI, hg, 
     hi, voli = 0}, Nnode = Length[coord]; If[Length[vol] == 0, 
      pvol = ConstantArray[vol, Nnode]; ]; If[Length[vol] == Nnode, 
      pvol = vol]; pwList = EmptyList[Length[nodeList]]; 
     hgList = EmptyList[Length[nodeList]]; 
     Do[NeiI = Nei[[i]]; coordi = ListMinus[coord[[NeiI]], coord[[i]]]; 
       voli = pvol[[NeiI]]; hi = Norm[coordi[[-1]]]; 
       {pw, khg} = Khgs[coordi, pfun, WeiF, voli, hgPen, hi]; 
       pw = PhuToPu[pw, iDfsum, iDffact, hi]; pwList[[kp]] = pw; 
       hgList[[kp++]] = khg; , {i, nodeList}]; {pwList, hgList}]
Khgs[v1_List, pfun_, wfun_, penalty_, h_] := 
   Module[{len = Length[pfun[v1[[1]], 1.]], num = Length[v1], k, p, trH = 0., 
     r1, p1, w1, pkp, sc, ssc, khg, pw0, pw}, 
    pkp = ConstantArray[0., {num, num}]; k = ConstantArray[0., {len, len}]; 
     p = ConstantArray[0., {num, len}]; 
     khg = ConstantArray[0., {num + 1, num + 1}]; 
     Do[r1 = Norm[v1[[i]]]; w1 = wfun[r1]; trH += w1*r1*r1; pkp[[i,i]] = w1; 
       p1 = pfun[v1[[i]], 1.1*h]; p[[i]] = w1*p1; 
       k += w1*p1 \[TensorProduct] p1, {i, num}]; 
     pw0 = Inverse[k] . Transpose[p]; pkp -= p . pw0; pkp *= penalty/trH; 
     sc = -Total[pkp, {1}]; ssc = -Total[sc]; khg[[1,1]] = ssc; 
     khg[[1,2 ;; num + 1]] = sc; khg[[2 ;; num + 1,1]] = sc; 
     khg[[2 ;; num + 1,2 ;; num + 1]] = pkp; 
     pw = ConstantArray[0., {len, num + 1}]; pw[[All,1]] = -Total[pw0, {2}]; 
     pw[[All,2 ;; num + 1]] = pw0; {pw, khg}]
 
Khgs[v1_List, pfun_, wfun_, vol_List, penalty_, h_] := 
   Module[{len = Length[pfun[v1[[1]], 1.]], num = Length[v1], k, p, trH = 0., 
     r1, p1, w1, pkp, sc, ssc, khg, pw0, pw, wl}, 
    pkp = ConstantArray[0., {num, num}]; k = ConstantArray[0., {len, len}]; 
     p = ConstantArray[0., {num, len}]; 
     khg = ConstantArray[0., {num + 1, num + 1}]; wl = Norm /@ v1; 
     Do[wl[[i]] = wfun[wl[[i]], 1.1*h]*vol[[i]], {i, num}]; 
     wl *= 1./Total[wl]; Do[r1 = Norm[v1[[i]]]; w1 = wl[[i]]; trH += w1; 
       pkp[[i,i]] = w1; p1 = pfun[v1[[i]], h]; p[[i]] = w1*p1; 
       k += w1*p1 \[TensorProduct] p1, {i, num}]; 
     pw0 = PseudoInverse[k] . Transpose[p]; pkp -= p . pw0; 
     pkp *= penalty/trH; sc = -Total[pkp, {1}]; ssc = -Total[sc]; 
     khg[[1,1]] = ssc; khg[[1,2 ;; num + 1]] = sc; khg[[2 ;; num + 1,1]] = 
      sc; khg[[2 ;; num + 1,2 ;; num + 1]] = pkp; 
     pw = ConstantArray[0., {len, num + 1}]; pw[[All,1]] = -Total[pw0, {2}]; 
     pw[[All,2 ;; num + 1]] = pw0; {pw, khg}]
PhuToPu[pw_List, iDfsum_List, iDffact_List, h_] := 
   Module[{p = {}, i}, p = pw; Do[p[[i]] *= 1./(h^iDfsum[[i]]/iDffact[[i]]), 
      {i, Length[pw]}]; p]
NOMResiPreValue2[NOMResi1_, ndim_Integer, pvol_List, nonvarList_List, 
    uvwN_List, PW_List, Nei2_List, DofNei_] := 
   Block[{uvwIlist, len = Length[pvol], udim = Length[uvwN[[1]]], nNei1 = 1, 
     Rvalue, Rv}, nNei1 = Length[Nei2[[1]]]; 
     uvwIlist = ConstantArray[0., {len, nNei1, udim}]; 
     Do[uvwIlist[[i]] = uvwN[[Prepend[Nei2[[i]], i]]], {i, len}]; 
     Rvalue = NOMResi1[ndim, pvol, nonvarList, uvwIlist, PW]; 
     Rv = ConstantArray[0., udim*len]; Do[Rv[[DofNei[[ni]]]] += Rvalue[[ni]], 
      {ni, len}]; Rv]
NOMResiPreValue3[NOMResi1_,ndim_Integer,pvol_List,nonvarList_List,uvwN_List,PW_List,KHG_List,weiNei_List,Nei2_List,DofNei_,phg_]:=Block[{uvwIlist,len=Length[pvol],udim=Length[uvwN[[1]]],nNei1=1,Rvalue,Rv,duvwi,hgContribution,abshg,hourglassForce},nNei1=Length[Nei2[[1]]];
uvwIlist=ConstantArray[0.,{len,nNei1,udim}];
Do[uvwIlist[[i]]=uvwN[[Prepend[Nei2[[i]],i]]],{i,len}];
duvwi=ConstantArray[0.,{len,udim,Length[PW[[1]]]+1}];
Do[duvwi[[i]]=Transpose[Prepend[PW[[i]].uvwIlist[[i]],uvwN[[i]]]],{i,len}];
(*xxx3=duvwi;*)
hgContribution=ConstantArray[0.,Dimensions[uvwIlist]];
Do[hgContribution[[i]]=KHG[[i]].uvwIlist[[i]],{i,len}];
abshg=AbsoluteHourglassValue[hgContribution,weiNei,pvol[[1]]];
(*xxx=abshg;*)
If[Length[phg]>1,Do[hgContribution[[i]]*=phg[[i]]*pvol[[1]],{i,len}],hgContribution*=phg*pvol[[1]];];
(*for variable phg, for the constraint of clamped boundary condition in thin plate;*)
(*hgC=hgContribution;*)
hourglassForce=Table[Flatten[hgContribution[[i]]],{i,len}];
Rvalue=NOMResi1[ndim,pvol,nonvarList,uvwIlist,duvwi,PW];
Rvalue+=hourglassForce;
Rv=ConstantArray[0.,udim len];
Do[Rv[[DofNei[[ni]]]]+=Rvalue[[ni]],{ni,len}];
Rv];
ParticleNeighDof[Nei_, udim_Integer] := 
   Table[Flatten[Transpose[Table[udim*(Prepend[Nei[[ni]], ni] - 1) + i, 
       {i, udim}]]], {ni, Length[Nei]}]
Plot2DField[vals_List, mesh_ElementMesh] := 
   Block[{vi, vm, color = "Rainbow", ntick = 9}, 
    {vm, vi} = EffectiveRange[vals, 0.025]; 
     Legended[Graphics[ElementMeshToGraphicsComplex[mesh, 
        VertexColors -> ColorData[color] /@ ((1./(vi - vm + 10.^(-20)))*
           (vals - vm))]], BarLegend[{color, {vm, vi}}, ntick]]]
 
Plot2DField[vals_List, mesh_ElementMesh, {vm_, vi_}] := 
   Block[{color = "Rainbow", ntick = 9}, 
    Legended[Graphics[ElementMeshToGraphicsComplex[mesh, 
       VertexColors -> ColorData[color] /@ ((1./(vi - vm + 10.^(-20)))*
          (vals - vm))]], BarLegend[{color, {vm, vi}}, ntick]]]
 
Plot2DField[vals_List, mesh_ElementMesh, uvw_List] := 
   Plot2DField[vals, DeformMesh[mesh, uvw]]
 
Plot2DField[vals_List, mesh_List] := 
   Block[{}, Plot2DField[vals, MyDelaunayMesh[mesh]]]
 
Plot2DField[vals_List, mesh_List, uvw_List] := 
   Block[{}, Plot2DField[vals, MyDelaunayMesh[mesh, uvw]]]
 
Plot2DField[vals_List, mesh_ElementMesh, ntick_Integer, tolRange_:0.001] := 
   Block[{vi, vm, color = "Rainbow", nInc}, 
    {vm, vi} = EffectiveRange[vals, tolRange]; 
     nInc = (vi - vm + 10.^(-50))/ntick; 
     Legended[Graphics[ElementMeshToGraphicsComplex[mesh, 
        VertexColors -> ColorData[color] /@ ((1./(vi - vm + 10.^(-50)))*
           (nInc*Round[(vals - vm)/nInc]))]], BarLegend[{color, {vm, vi}}, 
       Ticks -> Range[0, ntick]*nInc + vm, LegendMarkerSize -> 300]]]
EffectiveRange[a_List, tol_:0.025] := Block[{st, stl = 1, x1 = Length[a]*tol, 
     i1 = 1, i2 = 1, s1 = 0}, st = N[HistogramList[a]]; 
     stl = Length[st[[2]]]; Do[s1 += st[[2,i]]; If[s1 > x1, i1 = i; Break[]], 
      {i, stl}]; s1 = 0; Do[s1 += st[[2,i]]; If[s1 > x1, i2 = i; Break[]], 
      {i, stl, 1, -1}]; st[[1,{i1, i2 + 1}]]]
DeformMesh[mesh_ElementMesh, uvw_List] := 
   ToElementMesh["Coordinates" -> mesh["Coordinates"] + uvw, 
    "MeshElements" -> {mesh["MeshElements"][[1]]}]
ShowMesh[m_, withLable_:False] := If[withLable, 
    Show[m["Wireframe"["MeshElementIDStyle" -> Black]], 
     m["Wireframe"["MeshElement" -> "PointElements", "MeshElementIDStyle" -> 
        Red]]], Show[m["Wireframe"]]]
WeiNeiList[coord_, weiF_, Nei_] := Block[{len = Length[Nei], hi, 
     len2 = Length[Nei[[1]]], coordi, weiNei, weiNeiI, rij}, 
    weiNei = Table[0., {len}, {Length[Nei[[1]]] + 1}]; 
     Do[coordi = coord[[Nei[[i]]]]; weiNeiI = Table[0., {len2}]; 
       rij = coordi[[-1]] - coord[[i]]; hi = 1.1*Sqrt[rij . rij]; 
       Do[rij = coordi[[j]] - coord[[i]]; weiNeiI[[j]] = 
          weiF[Sqrt[rij . rij], hi]; , {j, len2}]; 
       weiNeiI *= 1/Total[weiNeiI]; weiNei[[i,2 ;; -1]] = weiNeiI, {i, len}]; 
     weiNei]
Wfr2[r_, h_] := 1./r^2
ListMatrixMultiply = CompiledFunction[{11, 12., 5852}, 
    {{_Real, 2}, {_Real, 2}}, {{3, 2, 0}, {3, 2, 1}, {3, 2, 2}}, 
    {{4, {2, 0, 0}}}, {0, 1, 0, 0, 3}, {{42, "Dot", 3, 2, 0, 3, 2, 1, 2, 0, 
      0, 3, 2, 2}, {1}}, Function[{a, b}, a . b, Listable], Evaluate]
pfun = CompiledFunction[{11, 12., 5468}, {{_Real, 1}, _Real}, 
    {{3, 1, 0}, {3, 0, 0}, {3, 1, 1}}, {{2, {2, 0, 1}}, {1, {2, 0, 2}}}, 
    {1, 3, 6, 0, 2}, {{40, 60, 3, 0, 0, 3, 0, 1}, {41, 259, 3, 0, 1, 3, 1, 0, 
      3, 1, 1}, {33, 1, 0}, {24, 0, 1, 0}, {32, 0, 0}, {2, 0, 3}, {49}, 
     {3, 1}, {37, 1, 2, 3, 1}, {37, 1, 1, 3, 2}, {40, 56, 3, 0, 2, 3, 0, 3}, 
     {16, 1, 2, 4}, {40, 56, 3, 0, 1, 3, 0, 5}, {34, 1, 5, 2, 1, 3, 4, 5, 3, 
      1}, {1}}, Function[{v1, h}, Block[{x, y}, {x, y} = v1/h; 
       {y, x, y^2, x*y, x^2}]], Evaluate]
AbsoluteHourglassValue = CompiledFunction[{11, 12., 7900}, 
    {{_Real, 2}, {_Real, 1}, _Real}, {{3, 2, 0}, {3, 1, 1}, {3, 0, 0}, 
     {3, 1, 4}}, {{0, {2, 0, 5}}, {4, {2, 0, 7}}, {1.*^-6, {3, 0, 5}}, 
     {1, {2, 0, 3}}, {7., {3, 0, 6}}, {0., {3, 0, 1}}}, {1, 8, 7, 0, 6}, 
    {{33, 1, 6}, {6, 5, 2}, {35, 6, 3, 4}, {6, 5, 4}, {3, 2}, 
     {36, 2, 1, 3, 4}, {4, 4, 6, -1}, {40, 60, 3, 0, 0, 3, 0, 2}, {33, 1, 2}, 
     {6, 5, 4}, {3, 15}, {38, 0, 0, 4, 1, 3}, {38, 1, 0, 4, 0, 3}, 
     {40, 38, 3, 0, 3, 3, 0, 4}, {27, 3, 6, 4, 5, 0}, {2, 0, 3}, {7, 1, 4}, 
     {3, 4}, {38, 1, 0, 4, 0, 4}, {40, 60, 3, 0, 4, 3, 0, 3}, {7, 3, 4}, 
     {41, 259, 3, 0, 4, 3, 1, 3, 3, 1, 5}, {42, "DotVV", 3, 1, 5, 3, 1, 5, 2, 
      0, 7, 3, 0, 4}, {40, 57, 3, 0, 4, 3, 0, 3}, {39, 4, 0, 4, 0, 3}, 
     {4, 4, 2, -14}, {1}}, Function[{hgContribution, weiNei, particleVolume}, 
     Block[{hg2 = Table[0., {Length[weiNei]}], hi, invp = 1/particleVolume}, 
      Do[hi = hgContribution[[i]]*If[Abs[weiNei[[i]]] < 1.*^-6, 0., 
            1/weiNei[[i]]]; hg2[[i]] = Sqrt[hi . hi], {i, Length[weiNei]}]; 
       hg2], Listable], Evaluate]
Resi1 = CompiledFunction[{11, 12., 7836}, {_Integer, _Real, {_Real, 1}, 
     {_Real, 2}, {_Real, 2}, {_Real, 2}}, {{2, 0, 0}, {3, 0, 0}, {3, 1, 0}, 
     {3, 2, 1}, {3, 2, 2}, {3, 2, 3}, {3, 1, 10}}, 
    {{0, {2, 0, 1}}, {11, {2, 0, 13}}, {{{4, 5, 6}}, {2, 2, 4}}, 
     {1., {3, 0, 25}}, {4, {2, 0, 7}}, {6, {2, 0, 9}}, {5, {2, 0, 8}}, 
     {2, {2, 0, 10}}, {-1, {2, 0, 26}}, {15, {2, 0, 2}}, {12, {2, 0, 14}}, 
     {7, {2, 0, 12}}, {1, {2, 0, 4}}, {3, {2, 0, 11}}, {14, {2, 0, 15}}, 
     {0., {3, 0, 1}}}, {1, 27, 26, 0, 14}, 
    {{6, 1, 3}, {42, "CopyTensor", 2, 2, 4, 2, 2, 5}, {6, 2, 6}, {6, 1, 16}, 
     {35, 6, 3, 6}, {6, 1, 5}, {3, 2}, {36, 16, 1, 3, 6}, {4, 5, 6, -1}, 
     {33, 5, 16}, {38, 2, 0, 4, 0, 7, 0, 4}, {38, 2, 0, 4, 0, 8, 0, 3}, 
     {38, 2, 0, 4, 0, 9, 0, 2}, {38, 0, 0, 4, 0, 7}, {38, 0, 0, 10, 0, 6}, 
     {38, 0, 0, 11, 0, 5}, {40, 56, 3, 0, 6, 3, 0, 8}, {39, 6, 0, 10, 0, 8}, 
     {38, 6, 0, 10, 0, 10}, {19, 10, 11}, {39, 6, 0, 11, 0, 11}, 
     {38, 6, 0, 11, 0, 10}, {10, 4, 12}, {13, 12, 10, 12}, 
     {39, 6, 0, 8, 0, 12}, {38, 6, 0, 8, 0, 10}, {40, 60, 3, 0, 10, 3, 0, 9}, 
     {39, 6, 0, 9, 0, 9}, {41, 263, 3, 0, 5, 2, 0, 11, 3, 0, 10}, 
     {39, 6, 0, 12, 0, 10}, {19, 6, 13}, {39, 6, 0, 13, 0, 13}, 
     {38, 6, 0, 13, 0, 15}, {10, 4, 14}, {13, 14, 15, 14}, 
     {39, 6, 0, 14, 0, 14}, {13, 4, 2, 15}, {39, 6, 0, 15, 0, 15}, 
     {38, 6, 0, 15, 0, 17}, {16, 6, 17, 19}, {39, 6, 0, 2, 0, 19}, 
     {38, 6, 0, 9, 0, 17}, {38, 6, 0, 12, 0, 16}, {16, 7, 17, 16, 4, 18}, 
     {10, 14, 17}, {40, 60, 3, 0, 17, 3, 0, 16}, {16, 18, 16, 18}, 
     {38, 6, 0, 9, 0, 16}, {38, 6, 0, 12, 0, 17}, {16, 7, 6, 16, 17, 2, 24}, 
     {10, 14, 16}, {40, 60, 3, 0, 16, 3, 0, 17}, {16, 24, 17, 24}, 
     {38, 6, 0, 9, 0, 17}, {38, 6, 0, 12, 0, 16}, {38, 6, 0, 14, 0, 21}, 
     {16, 21, 4, 21}, {38, 6, 0, 2, 0, 22}, {13, 21, 22, 21}, 
     {16, 7, 17, 16, 21, 22}, {10, 14, 17}, {40, 60, 3, 0, 17, 3, 0, 16}, 
     {16, 22, 16, 22}, {13, 18, 24, 22, 18}, {10, 10, 24}, 
     {40, 60, 3, 0, 24, 3, 0, 22}, {16, 18, 22, 18}, {38, 6, 0, 14, 0, 22}, 
     {38, 6, 0, 9, 0, 24}, {38, 6, 0, 12, 0, 16}, {16, 7, 22, 24, 16, 3, 17}, 
     {10, 9, 22}, {40, 60, 3, 0, 22, 3, 0, 24}, {16, 17, 24, 17}, 
     {38, 6, 0, 9, 0, 24}, {38, 6, 0, 12, 0, 22}, {16, 7, 6, 24, 22, 4, 16}, 
     {10, 14, 24}, {40, 60, 3, 0, 24, 3, 0, 22}, {16, 16, 22, 16}, 
     {38, 6, 0, 9, 0, 22}, {38, 6, 0, 12, 0, 24}, {16, 7, 22, 24, 2, 21}, 
     {10, 14, 22}, {40, 60, 3, 0, 22, 3, 0, 24}, {16, 21, 24, 21}, 
     {38, 6, 0, 9, 0, 24}, {38, 6, 0, 12, 0, 22}, {38, 6, 0, 14, 0, 20}, 
     {16, 20, 2, 20}, {38, 6, 0, 2, 0, 23}, {13, 20, 23, 20}, 
     {16, 7, 24, 22, 20, 23}, {10, 14, 24}, {40, 60, 3, 0, 24, 3, 0, 22}, 
     {16, 23, 22, 23}, {13, 16, 21, 23, 16}, {10, 10, 21}, 
     {40, 60, 3, 0, 21, 3, 0, 23}, {16, 16, 23, 16}, 
     {34, 1, 3, 18, 17, 16, 3, 7}, {38, 3, 0, 4, 1, 8}, {33, 8, 5}, 
     {33, 5, 18}, {6, 1, 19}, {3, 14}, {38, 5, 0, 19, 1, 10}, {33, 10, 20}, 
     {6, 1, 22}, {3, 9}, {37, 10, 22, 2, 23}, {24, 23, 1, 0}, {2, 0, 3}, 
     {3, 6}, {3, 1}, {6, 3, 24}, {12, 24, 4, 25}, {6, 25, 3}, 
     {4, 22, 20, -8}, {4, 19, 18, -13}, {6, 3, 17}, {15, 16, 5, 24}, 
     {6, 1, 23}, {35, 17, 24, 3, 8}, {6, 1, 20}, {3, 5}, {6, 1, 22}, {3, 2}, 
     {36, 23, 1, 3, 8}, {4, 22, 24, -1}, {4, 20, 17, -4}, {6, 4, 6}, 
     {6, 16, 18}, {6, 1, 19}, {3, 26}, {38, 5, 0, 19, 1, 11}, {33, 11, 23}, 
     {6, 1, 22}, {3, 21}, {37, 11, 22, 2, 20}, {24, 20, 1, 0}, {2, 0, 3}, 
     {3, 18}, {3, 1}, {24, 20, 4, 0}, {32, 0, 0}, {2, 0, 9}, {6, 6, 17}, 
     {12, 17, 4, 24}, {6, 24, 6}, {34, 1, 4, 19, 26, 16, 1, 2, 12}, 
     {12, 20, 26, 25}, {38, 3, 0, 25, 1, 13}, {39, 8, 0, 17, 3, 12, 1, 13}, 
     {3, 5}, {6, 6, 17}, {12, 17, 4, 25}, {6, 25, 6}, 
     {39, 8, 0, 17, 0, 19, 0, 25}, {4, 22, 23, -20}, {4, 19, 18, -25}, 
     {41, 259, 3, 0, 0, 3, 1, 7, 3, 1, 9}, {42, "Dot", 3, 1, 9, 3, 2, 8, 2, 
      0, 7, 3, 1, 10}, {1}}, Function[{ndim, volSet, nonvars, uvwi, duvwi, 
      pw}, Block[{udim, tnablaU, NeiI2, Ri, k, ulen = 0, nNode, 
       tidf2 = {{4, 5, 6}}, $v = Table[0., {15}], El$$, nu$$, thickness$$, 
       u02, u11, u20}, udim = Length[tidf2]; u02 = duvwi[[1,4]]; 
       u11 = duvwi[[1,5]]; u20 = duvwi[[1,6]]; El$$ = nonvars[[1]]; 
       nu$$ = nonvars[[2]]; thickness$$ = nonvars[[3]]; $v[[2]] = nu$$^2; 
       $v[[3]] = -$v[[2]]; $v[[5]] = 1 + $v[[3]]; $v[[6]] = $v[[5]]^(-1); 
       $v[[7]] = thickness$$^3; $v[[11]] = -nu$$; $v[[12]] = 1 + $v[[11]]; 
       $v[[14]] = u02 + u20; $v[[15]] = nu$$*$v[[14]]; 
       Ri = {((El$$*$v[[6]]*$v[[7]]*u02)/12 + (El$$*nu$$*$v[[6]]*$v[[7]]*u20)/
            12 + (El$$*$v[[6]]*$v[[7]]*($v[[12]]*u02 + $v[[15]]))/12)/2, 
         (El$$*$v[[12]]*$v[[6]]*$v[[7]]*u11)/6, 
         ((El$$*nu$$*$v[[6]]*$v[[7]]*u02)/12 + (El$$*$v[[6]]*$v[[7]]*u20)/
            12 + (El$$*$v[[6]]*$v[[7]]*($v[[12]]*u20 + $v[[15]]))/12)/2}; 
       nNode = Length[pw[[1]]]; Do[Do[If[j == 0, Break[]]; ulen++, 
         {j, tidf2[[i]]}], {i, Length[tidf2]}]; 
       tnablaU = Table[0., {ulen}, {udim*nNode}]; k = 1; 
       Do[Do[If[j == 0, Break[]]; If[j != 1, tnablaU[[k++,i ;; -1 ;; udim]] = 
            pw[[j - 1]], tnablaU[[k++,i]] = 1.], {j, tidf2[[i]]}], 
        {i, udim}]; (volSet*Ri) . tnablaU], Listable], Evaluate]
Resi2 = CompiledFunction[{11, 12., 7836}, {_Integer, _Real, {_Real, 1}, 
     {_Real, 2}, {_Real, 2}, {_Real, 2}}, {{2, 0, 0}, {3, 0, 0}, {3, 1, 0}, 
     {3, 2, 1}, {3, 2, 2}, {3, 2, 3}, {3, 1, 10}}, 
    {{0, {2, 0, 1}}, {{{4, 5, 6}}, {2, 2, 4}}, {1., {3, 0, 23}}, 
     {9, {2, 0, 13}}, {4, {2, 0, 7}}, {8, {2, 0, 12}}, {13, {2, 0, 2}}, 
     {6, {2, 0, 9}}, {5, {2, 0, 8}}, {2, {2, 0, 10}}, {-1, {2, 0, 26}}, 
     {12, {2, 0, 15}}, {7, {2, 0, 14}}, {1, {2, 0, 4}}, {3, {2, 0, 11}}, 
     {0., {3, 0, 1}}}, {1, 27, 24, 0, 14}, 
    {{6, 1, 3}, {42, "CopyTensor", 2, 2, 4, 2, 2, 5}, {6, 2, 6}, {6, 1, 16}, 
     {35, 6, 3, 6}, {6, 1, 5}, {3, 2}, {36, 16, 1, 3, 6}, {4, 5, 6, -1}, 
     {33, 5, 16}, {38, 2, 0, 4, 0, 7, 0, 4}, {38, 2, 0, 4, 0, 8, 0, 3}, 
     {38, 2, 0, 4, 0, 9, 0, 2}, {38, 0, 0, 4, 0, 7}, {38, 0, 0, 10, 0, 6}, 
     {38, 0, 0, 11, 0, 5}, {19, 6, 8}, {39, 6, 0, 12, 0, 8}, 
     {38, 6, 0, 12, 0, 10}, {10, 4, 9}, {13, 9, 10, 9}, {39, 6, 0, 13, 0, 9}, 
     {40, 56, 3, 0, 6, 3, 0, 10}, {39, 6, 0, 10, 0, 10}, 
     {38, 6, 0, 10, 0, 12}, {19, 12, 13}, {39, 6, 0, 11, 0, 13}, 
     {38, 6, 0, 11, 0, 12}, {10, 4, 14}, {13, 14, 12, 14}, 
     {39, 6, 0, 8, 0, 14}, {38, 6, 0, 8, 0, 12}, {40, 60, 3, 0, 12, 3, 0, 
      11}, {39, 6, 0, 9, 0, 11}, {41, 263, 3, 0, 5, 2, 0, 11, 3, 0, 12}, 
     {39, 6, 0, 14, 0, 12}, {13, 4, 2, 15}, {39, 6, 0, 15, 0, 15}, 
     {38, 6, 0, 15, 0, 17}, {16, 6, 17, 20}, {39, 6, 0, 2, 0, 20}, 
     {38, 6, 0, 9, 0, 17}, {38, 6, 0, 14, 0, 21}, {38, 6, 0, 13, 0, 22}, 
     {16, 22, 4, 22}, {38, 6, 0, 2, 0, 18}, {13, 22, 18, 22}, 
     {16, 7, 17, 21, 22, 18}, {10, 15, 17}, {40, 60, 3, 0, 17, 3, 0, 21}, 
     {16, 18, 21, 18}, {38, 6, 0, 13, 0, 21}, {38, 6, 0, 9, 0, 17}, 
     {38, 6, 0, 14, 0, 22}, {16, 7, 21, 17, 22, 3, 16}, {10, 9, 21}, 
     {40, 60, 3, 0, 21, 3, 0, 17}, {16, 16, 17, 16}, {38, 6, 0, 9, 0, 17}, 
     {38, 6, 0, 14, 0, 21}, {38, 6, 0, 13, 0, 22}, {16, 22, 2, 22}, 
     {38, 6, 0, 2, 0, 19}, {13, 22, 19, 22}, {16, 7, 17, 21, 22, 19}, 
     {10, 15, 17}, {40, 60, 3, 0, 17, 3, 0, 21}, {16, 19, 21, 19}, 
     {34, 1, 3, 18, 16, 19, 3, 7}, {38, 3, 0, 4, 1, 8}, {33, 8, 5}, 
     {33, 5, 18}, {6, 1, 19}, {3, 14}, {38, 5, 0, 19, 1, 10}, {33, 10, 20}, 
     {6, 1, 22}, {3, 9}, {37, 10, 22, 2, 23}, {24, 23, 1, 0}, {2, 0, 3}, 
     {3, 6}, {3, 1}, {6, 3, 24}, {12, 24, 4, 25}, {6, 25, 3}, 
     {4, 22, 20, -8}, {4, 19, 18, -13}, {6, 3, 17}, {15, 16, 5, 24}, 
     {6, 1, 23}, {35, 17, 24, 3, 8}, {6, 1, 20}, {3, 5}, {6, 1, 22}, {3, 2}, 
     {36, 23, 1, 3, 8}, {4, 22, 24, -1}, {4, 20, 17, -4}, {6, 4, 6}, 
     {6, 16, 18}, {6, 1, 19}, {3, 26}, {38, 5, 0, 19, 1, 11}, {33, 11, 23}, 
     {6, 1, 22}, {3, 21}, {37, 11, 22, 2, 20}, {24, 20, 1, 0}, {2, 0, 3}, 
     {3, 18}, {3, 1}, {24, 20, 4, 0}, {32, 0, 0}, {2, 0, 9}, {6, 6, 17}, 
     {12, 17, 4, 24}, {6, 24, 6}, {34, 1, 4, 19, 26, 16, 1, 2, 12}, 
     {12, 20, 26, 25}, {38, 3, 0, 25, 1, 13}, {39, 8, 0, 17, 3, 12, 1, 13}, 
     {3, 5}, {6, 6, 17}, {12, 17, 4, 25}, {6, 25, 6}, 
     {39, 8, 0, 17, 0, 19, 0, 23}, {4, 22, 23, -20}, {4, 19, 18, -25}, 
     {41, 259, 3, 0, 0, 3, 1, 7, 3, 1, 9}, {42, "Dot", 3, 1, 9, 3, 2, 8, 2, 
      0, 7, 3, 1, 10}, {1}}, Function[{ndim, volSet, nonvars, uvwi, duvwi, 
      pw}, Block[{udim, tnablaU, NeiI2, Ri, k, ulen = 0, nNode, 
       tidf2 = {{4, 5, 6}}, $v = Table[0., {13}], El$$, nu$$, thickness$$, 
       u02, u11, u20}, udim = Length[tidf2]; u02 = duvwi[[1,4]]; 
       u11 = duvwi[[1,5]]; u20 = duvwi[[1,6]]; El$$ = nonvars[[1]]; 
       nu$$ = nonvars[[2]]; thickness$$ = nonvars[[3]]; $v[[8]] = -nu$$; 
       $v[[9]] = 1 + $v[[8]]; $v[[2]] = nu$$^2; $v[[3]] = -$v[[2]]; 
       $v[[5]] = 1 + $v[[3]]; $v[[6]] = $v[[5]]^(-1); 
       $v[[7]] = thickness$$^3; $v[[12]] = u02 + u20; 
       $v[[13]] = nu$$*$v[[12]]; 
       Ri = {(El$$*$v[[6]]*$v[[7]]*($v[[9]]*u02 + $v[[13]]))/12, 
         (El$$*$v[[9]]*$v[[6]]*$v[[7]]*u11)/6, 
         (El$$*$v[[6]]*$v[[7]]*($v[[9]]*u20 + $v[[13]]))/12}; 
       nNode = Length[pw[[1]]]; Do[Do[If[j == 0, Break[]]; ulen++, 
         {j, tidf2[[i]]}], {i, Length[tidf2]}]; 
       tnablaU = Table[0., {ulen}, {udim*nNode}]; k = 1; 
       Do[Do[If[j == 0, Break[]]; If[j != 1, tnablaU[[k++,i ;; -1 ;; udim]] = 
            pw[[j - 1]], tnablaU[[k++,i]] = 1.], {j, tidf2[[i]]}], 
        {i, udim}]; (volSet*Ri) . tnablaU], Listable], Evaluate]
        
(*Plot the deflection (u) of the points in 2D (coord)*)
MyPlotPoint[coord_List,u_List,plotRange_:Automatic]:=Module[{Nnode=Length[coord],ndim=Length[coord[[1]]],xyz={},uv={"u","v","w"},r,r2},
If[Length[coord]>Length[u],Print["Error,u is less than coord"];Abort[]];xyz=ConstantArray[0.,{Nnode,ndim+1}];xyz[[All,1;;ndim]]=coord;
If[Length[u]==Nnode,xyz[[All,ndim+1]]=u;CellPrint@ExpressionCell[#,"Output"]&/@{ListPointPlot3D[xyz,PlotLegends->Automatic, Axes->True,AxesLabel->{"x","y"},PlotLabel->"u",PlotRange->plotRange, ColorFunction->"SouthwestColors"]},
Do[xyz[[All,ndim+1]]=u[[iu;;-1;;ndim]];
CellPrint@ExpressionCell[#,"Output"]&/@{ListPointPlot3D[xyz,PlotLegends->Automatic,Axes->True,AxesLabel->{"x","y"},PlotLabel->uv[[iu]], ColorFunction->"Rainbow"]},{iu,ndim}]]];
MyPlotPoint[coord_List,u_List,pset_List,plotRange_:Automatic]:=Module[{Nnode=Length[coord],ndim=Length[coord[[1]]],xyz={},uv={"u","v","w"},r,r2},
If[Length[coord]>Length[u],Print["Error,u is less than coord"];Abort[]];xyz=ConstantArray[0.,{Nnode,ndim+1}];xyz[[All,1;;ndim]]=coord;
If[Length[u]==Nnode,xyz[[All,ndim+1]]=u;CellPrint@ExpressionCell[#,"Output"]&/@{ListPointPlot3D[xyz[[pset,All]],PlotLegends->Automatic, Axes->True,AxesLabel->{"x","y"},PlotLabel->"u",PlotRange->plotRange, ColorFunction->"SouthwestColors"]},
Do[xyz[[All,ndim+1]]=u[[iu;;-1;;ndim]];
CellPrint@ExpressionCell[#,"Output"]&/@{ListPointPlot3D[xyz[[pset]],PlotLegends->Automatic,Axes->True,AxesLabel->{"x","y"},PlotLabel->uv[[iu]], ColorFunction->"Rainbow"]},{iu,ndim}]]];

ClearAll[SaveToCell]
SaveToCell::usage="SaveToCell[variable] creates an input cell that reassigns the current value of variable.\n"<>"SaveToCell[variables, display] shows 'display' on the right-hand-side of the assignment.";
SetAttributes[SaveToCell,HoldFirst]
SaveToCell[var_,name:Except[_?OptionQ]:"data",opt:OptionsPattern[]]:=With[{data=Compress[var],panel=ToBoxes@Tooltip[Panel[name,FrameMargins->Small],DateString[]]},CellPrint@Cell[BoxData@RowBox[{MakeBoxes[var],"=",InterpretationBox[panel,Uncompress[data]],";"}],"Input",(*prevent deletion by Cell>Delete All Output:*)GeneratedCell->False,(*CellLabel is special:last occrrence takes precedence,so it comes before opt:*)CellLabel->"(saved)",opt,CellLabelAutoDelete->False]];
