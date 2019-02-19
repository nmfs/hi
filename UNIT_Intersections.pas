unit UNIT_Intersections;

interface

uses

    math,
    sysUtils,
    Unit_Utils,
    Unit_Clusters,
    Unit_Parameterizations,
    Mesh_Front_1D
    ;

type



  Func1D = function(const x: double): double;

  FuncAreWithinTol = function(const t1, t2, tol: double): boolean of Object;


  TObjectType = (otPoint, otCurve, otSurface, otSolid);


  TObject3D = class;


  TObjectArray = array of TObject3D;


  TProtoObject3D = record

    Obj: TObject3D;
    ParaPoint: array of double;
    Inward: array of boolean;

    function INWARD_INDEX(const inward: boolean): integer;
  end;


  TProtoObjectArray = array of TProtoObject3D;



  TObject3D = class

    ObjType: TObjectType;
    Param: TParam;
    Boundary: TObjectArray;
    ParaPoint: array of double;
    Inward: array of boolean;

    Contain: TProtoObjectArray;
    ContainedBy: TObjectArray;

    Cluster: array of TCluster;

    net: TNet;


    procedure GENERATE_BASE_CLUSTERS(const safeStep: double); Virtual; Abstract;

    procedure ADD_CONTAINED_BY(const obj: TObject3d);

    function INWARD_INDEX(const inward: boolean): integer;

  end;




  TPoint = class(TObject3D)

    constructor Create;

    procedure GENERATE_BASE_CLUSTERS(const safeStep: double); Override;

  end;


  TCurve = class(TObject3D)

    startPoint, endPoint: TParametricCurvePoint;
    isLoop: boolean;

    constructor Create;

    procedure GENERATE_BASE_CLUSTERS(const safeStep: double); Override;

    function IS_INSIDE(const t: double; const tol: double): boolean;

  end;


  TSurface = class(TObject3D)

    constructor Create;

  end;


  TSolid = class(TObject3D)

    constructor Create;

  end;


  TGroups = record

    group: array of TObjectArray;

    procedure ASSOCIATION(const obj: TObjectArray; const tolL: double);
    procedure ADD_TO_GROUP(const obj1, obj2: TObject3d);
    function GROUP_INDEX(const obj: TObject3d): integer;
    function ARE_ASSOCIATED(const obj1, obj2: TObject3d): boolean;
    procedure REDUCE(const obj: TObjectArray);

  end;



  TGmshSegment = record

    name: integer;
    nodei: TVec3d;
    nodej: TVec3d;
    namei: integer;
    namej: integer;
  end;

  TGmshState = record

    first: integer;
    last: integer;
    firstP: integer;
    lastP: integer;
  end;

  TGmshUnique = record

    name: integer;
    color: double;
    rand: TVec3d;
  end;

  TGmshFile = record

    points: array of TGmshSegment;
    segments: array of TGmshSegment;
    state: array of TGmshState;
    unique: array of TGmshUnique;

    procedure ADD_EXPLICIT(const obj: TObjectArray);
    procedure ADD_STATE(const obj: TObjectArray);
    procedure SAVE_TO_FILE;

    function GET_INDEX(const name: integer; const isPoint: boolean): integer;
    function GET_COLOR(const name: integer; const isPoint: boolean): double;
    function GET_RAND(const name: integer; const isPoint: boolean): TVec3d;
    procedure FIND_UNIQUE;
  end;


  //============================================================================
  //FUNCTIONS
  //============================================================================

  function INTERSECT(const obj1, obj2: TObject3D;
                     const tol: double = 1e-3;
                     const tolL: double = 1e-3;
                     const tolA: double = 1e-3
                     ): TObjectArray; Overload;


  function RECONCILE(const obj: TObjectArray;
                     const tol: double = 1e-3;
                     const tolL: double = 1e-3;
                     const tolA: double = 1e-3
                     ): TObjectArray; Overload;



implementation

//============================================================================

procedure TObject3d.ADD_CONTAINED_BY(const obj: TObject3d);
var

  k: integer;
  bol: boolean;
begin

  bol:= false;
  for k:= 0 to length(ContainedBy)-1 do begin

    if obj = ContainedBy[k] then begin

      bol:= true;
    end;
  end;
  if bol = false then begin

    setlength(ContainedBy, k+1);
    ContainedBy[k]:= obj;
  end;
end;
//============================================================================

function TObject3d.INWARD_INDEX(const inward: boolean): integer;
begin

  if self.Inward[0] = inward then begin

    result:= 0;
  end else begin

    result:= 1;
  end;
end;
//============================================================================
function TProtoObject3D.INWARD_INDEX(const inward: boolean): integer;
begin

  if self.Inward[0] = inward then begin

    result:= 0;
  end else begin

    result:= 1;
  end;
end;
//============================================================================

procedure TGroups.ASSOCIATION(const obj: TObjectArray; const tolL: double);
var

  i, j, k: integer;
  rt: double;
begin

  for i:= 0 to length(obj)-1 do begin

    if (obj[i].ObjType = otCurve) and (obj[i].Param <> nil) then begin

      rt:= NORM(obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(0).rt);
    end;

    for j:= 0 to length(obj[i].Contain)-1 do begin

      if obj[i].Contain[j].Obj.ObjType = otCurve then begin

        for k:= j+1 to length(obj[i].Contain)-1 do begin

          if obj[i].Contain[k].Obj.ObjType = otCurve then begin

            if abs(obj[i].Contain[j].ParaPoint[0]-
                   obj[i].Contain[k].ParaPoint[0])*rt < tolL then begin

              ADD_TO_GROUP(obj[i].Contain[j].Obj.Boundary[0],
                           obj[i].Contain[k].Obj.Boundary[0]);
            end;
            if abs(obj[i].Contain[j].ParaPoint[0]-
                   obj[i].Contain[k].ParaPoint[1])*rt < tolL then begin

              ADD_TO_GROUP(obj[i].Contain[j].Obj.Boundary[0],
                           obj[i].Contain[k].Obj.Boundary[1]);
            end;
            if abs(obj[i].Contain[j].ParaPoint[1]-
                   obj[i].Contain[k].ParaPoint[0])*rt < tolL then begin

              ADD_TO_GROUP(obj[i].Contain[j].Obj.Boundary[1],
                           obj[i].Contain[k].Obj.Boundary[0]);
            end;
            if abs(obj[i].Contain[j].ParaPoint[1]-
                   obj[i].Contain[k].ParaPoint[1])*rt < tolL then begin

              ADD_TO_GROUP(obj[i].Contain[j].Obj.Boundary[1],
                           obj[i].Contain[k].Obj.Boundary[1]);
            end;
          end;
        end;
        if abs(obj[i].Contain[j].ParaPoint[0]-
               obj[i].Contain[j].ParaPoint[1])*rt < tolL then begin

          ADD_TO_GROUP(obj[i].Contain[j].Obj.Boundary[0],
                       obj[i].Contain[j].Obj.Boundary[1]);
        end;
      end;
    end;
  end;
end;
//============================================================================

procedure TGroups.REDUCE(const obj: TObjectArray);
var

  i, j, k: integer;
  rt: double;
begin

  for i:= 0 to length(obj)-1 do begin

    if (obj[i].ObjType = otCurve) and (obj[i].Param <> nil) then begin

      rt:= NORM(obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(0).rt);
    end;

    for j:= 0 to length(obj[i].Contain)-1 do begin

      if obj[i].Contain[j].Obj.ObjType = otCurve then begin

        k:= GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[0]);
        if k <> -1 then begin

          obj[i].Contain[j].Obj.Boundary[0]:= group[k][0];
        end;
      end;
    end;
  end;
end;
//==============================================================================

procedure TGroups.ADD_TO_GROUP(const obj1, obj2: TObject3d);
var

  k, kk, kkk, kkkk, i: integer;
begin

  k:= GROUP_INDEX(obj1);
  kk:= GROUP_INDEX(obj2);
  if (k = -1) and (kk = -1) then begin

    k:= length(group);
    setlength(group, k+1);
    setlength(group[k], 2);
    group[k][0]:= obj1;
    group[k][1]:= obj2;
  end else begin

    if k = -1 then begin

      kkk:= length(group[kk]);
      setlength(group[kk], kkk+1);
      group[kk][kkk]:= obj1;
    end;
    if kk = -1 then begin

      kkk:= length(group[k]);
      setlength(group[k], kkk+1);
      group[k][kkk]:= obj2;
    end;
    if (k <> -1) and (kk <> -1) and (k <> kk) then begin

      // merge groups
      kkk:= length(group[k]);
      kkkk:= length(group[kk]);
      setlength(group[k], kkk+kkkk);
      for i:= 0 to kkkk-1 do begin

        group[k][kkk+i]:= group[kk][i];
      end;
      setlength(group[kk], 0);
      kkk:= length(group);
      for i:= kk to kkk-2 do begin

        group[i]:= group[i+1];
      end;
      setlength(group, kkk-1);
    end;
  end;
end;
//==============================================================================

function TGroups.GROUP_INDEX(const obj: TObject3d): integer;
var

  k, kk: integer;
begin

  result:= -1;
  for k:= 0 to length(group)-1 do begin

    for kk:= 0 to length(group[k])-1 do begin

      if obj = group[k][kk] then begin

        result:= k;
      end;
    end;
  end;
end;
//==============================================================================
function TGroups.ARE_ASSOCIATED(const obj1, obj2: TObject3d): boolean;
var

  i, j: integer;
begin

  i:= GROUP_INDEX(obj1);
  j:= GROUP_INDEX(obj2);

  if (i = -1) and (j = -1) then begin

    result:= false;
  end else if i = j then begin

    result:= true;
  end else begin

    result:= false;
  end;
end;
//==============================================================================
function TGmshFile.GET_INDEX(const name: integer; const isPoint: boolean): integer;
var

  k: integer;
  value: double;
begin

  result:= -1;
  for k:= 0 to length(unique)-1 do begin

    if name = unique[k].name then begin

      result:= k;
    end;
  end;
  if result = -1 then begin

    k:= length(unique);
    setlength(unique, k+1);
    unique[k].name:= name;
    unique[k].color:= random;
    unique[k].rand[0]:= random;
    unique[k].rand[1]:= random;
    unique[k].rand[2]:= random;
    value:= 1e-3/NORM(unique[k].rand);
    unique[k].rand[0]:= value*unique[k].rand[0];
    unique[k].rand[1]:= value*unique[k].rand[1];
    unique[k].rand[2]:= value*unique[k].rand[2];
    if length(state) > 0 then begin

      if k = 0+state[0].first+state[0].firstP then unique[k].color:= 0.2;
      if k = 1+state[0].first+state[0].firstP then unique[k].color:= 0.6;
      if k = 2+state[0].first+state[0].firstP then unique[k].color:= 0.9;
      if k = 3+state[0].first+state[0].firstP then unique[k].color:= 0.4;
      if k = 4+state[0].first+state[0].firstP then unique[k].color:= 0.7;
      if k = 5+state[0].first+state[0].firstP then unique[k].color:= 0;
      if k = 6+state[0].first+state[0].firstP then unique[k].color:= 1;
      if k = 7+state[0].first+state[0].firstP then unique[k].color:= 0.45;
      if k = 8+state[0].first+state[0].firstP then unique[k].color:= 0.8;
      if k = 9+state[0].first+state[0].firstP then unique[k].color:= 0.3;
    end;
    result:= k;
  end;
end;
//==============================================================================

function TGmshFile.GET_COLOR(const name: integer; const isPoint: boolean): double;
begin

  result:= unique[GET_INDEX(name, isPoint)].color;
end;
//==============================================================================

function TGmshFile.GET_RAND(const name: integer; const isPoint: boolean): TVec3d;
begin

  result:= unique[GET_INDEX(name, isPoint)].rand;
end;
//==============================================================================
procedure TGmshFile.FIND_UNIQUE;
var

  k, kk, index: integer;
begin
   {
  setlength(unique, 0);
  for k:= 0 to length(points)-1 do begin

    index:=


      if name = points[k].name then begin

        result:= k;
      end;
  end;     }
end;
//==============================================================================

procedure TGmshFile.ADD_EXPLICIT(const obj: TObjectArray);
var

  i, j, k, kk, kkk, divisions: integer;
  incr: double;
begin

  setlength(segments, 0);
  setlength(points, 0);
  setlength(state, 0);
  setlength(unique, 0);
  divisions:= 10; // minimum 1 division
  for i:= 0 to length(obj)-1 do begin

    if obj[i].Param <> nil then begin

      if obj[i].ObjType = otPoint then begin

        kk:= length(points);
        setlength(points, kk+1);
        points[kk].name:= Integer(Pointer(obj[i]));
        points[kk].nodei:= obj[i].Param.EVALUATE;
      end else begin

        if obj[i].Param.ParamType = ptStraightLine then begin

          kk:= length(segments);
          setlength(segments, kk+1);
          segments[kk].name:= Integer(Pointer(obj[i]));
          segments[kk].namei:= Integer(Pointer(obj[i].Boundary[0]));
          segments[kk].namej:= Integer(Pointer(obj[i].Boundary[1]));
          //if length(obj[i].ParaPoint) > 1 then begin

            segments[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].ParaPoint[0]).r;
            segments[kk].nodej:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].ParaPoint[1]).r;
          //end;
        end else begin

          kkk:= length(obj[i].Cluster)-2;
          kk:= length(segments);
          setlength(segments, kk+kkk*divisions);
          for k:= 0 to kkk-1 do begin

            incr:= (obj[i].net[k+1].p[0]-obj[i].net[k].p[0])/divisions;
            for j:= 0 to divisions-1 do begin

              segments[kk+k*divisions+j].name:= Integer(Pointer(obj[i]));
              segments[kk+k*divisions+j].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].net[k].p[0]+j*incr).r;
              segments[kk+k*divisions+j].nodej:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].net[k].p[0]+(j+1)*incr).r;
            end;
          end;
          if obj[i].net[0].neighbor[0] <> -1 then begin

            // is loop
            kk:= length(segments);
            setlength(segments, kk+divisions);
            incr:= (obj[i].net[1].p[0]-obj[i].net[0].p[0])/divisions;
            for j:= 0 to divisions-1 do begin

              segments[kk+j].name:= Integer(Pointer(obj[i]));
              segments[kk+j].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].net[kkk].p[0]+j*incr).r;
              segments[kk+j].nodej:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].net[kkk].p[0]+(j+1)*incr).r;
            end;
          end;
        end;
      end;
    end;
  end;
end;
//==============================================================================

procedure TGmshFile.ADD_STATE(const obj: TObjectArray);
var

  k, kk, kkk, i, j: integer;
  obj2: TObject3d;
  color: array of TVec3d;
  t, t0, t1, incr: double;
  bol: boolean;
begin

  kkk:= length(state);
  setlength(state, kkk+1);
  state[kkk].first:= length(segments);
  state[kkk].firstP:= length(points);
  for i:= 0 to length(obj)-1 do begin

    if obj[i].Param <> nil then begin

      for j:= 0 to length(obj[i].Contain)-1 do begin

        obj2:= obj[i].Contain[j].Obj;
        if obj2 <> nil then begin

          if obj2.ObjType = otPoint then begin

            kk:= length(points);
            setlength(points, kk+1);
            points[kk].name:= Integer(Pointer(obj2));
            if obj[i].ObjType = otPoint then begin

              points[kk].nodei:= obj[i].Param.EVALUATE;
            end else begin

              points[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].Contain[j].ParaPoint[0]).r;
            end;
          end else begin

            t0:= obj[i].Contain[j].ParaPoint[0];
            t1:= obj[i].Contain[j].ParaPoint[1];
            if obj[i].Param.ParamType = ptStraightLine then begin

              kk:= length(segments);
              setlength(segments, kk+1);
              segments[kk].name:= Integer(Pointer(obj2));
              segments[kk].namei:= Integer(Pointer(obj2.Boundary[0]));
              segments[kk].namej:= Integer(Pointer(obj2.Boundary[1]));
              segments[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t0).r;
              segments[kk].nodej:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t1).r;
              kk:= length(points);
              setlength(points, kk+1);
              points[kk].name:= Integer(Pointer(obj2.Boundary[0]));
              points[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t0).r;
              kk:= length(points);
              setlength(points, kk+1);
              points[kk].name:= Integer(Pointer(obj2.Boundary[1]));
              points[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t1).r;
            end else begin

              incr:= 0.1*obj[i].Param.SAFE_STEP/NORM(obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t0).rt);
              if t0 < t1 then begin

                t:= t0;
              end else begin

                t:= t1;
              end;
              bol:= true;
              while bol do begin

                kk:= length(segments);
                setlength(segments, kk+1);
                segments[kk].name:= Integer(Pointer(obj2));
                segments[kk].namei:= -1;
                segments[kk].namej:= -1;
                segments[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t).r;
                t:= t+incr;
                if t0 < t1 then begin

                  if t >= t1 then begin

                    t:= t1;
                    bol:= false;
                  end;
                end else begin

                  if t >= t0 then begin

                    t:= t0;
                    bol:= false;
                  end;
                end;
                segments[kk].nodej:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t).r;
                kk:= length(points);
                setlength(points, kk+1);
                points[kk].name:= Integer(Pointer(obj2.Boundary[0]));
                points[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].Contain[j].ParaPoint[0]).r;
                kk:= length(points);
                setlength(points, kk+1);
                points[kk].name:= Integer(Pointer(obj2.Boundary[1]));
                points[kk].nodei:= obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(obj[i].Contain[j].ParaPoint[1]).r;
              end;
            end;
          end;
        end;
      end;
    end;
  end;
  state[kkk].last:= length(segments)-1;
  state[kkk].lastP:= length(points)-1;
end;
//==============================================================================

procedure TGmshFile.SAVE_TO_FILE;
var

  k, kk: integer;
  myFile: TextFile;
  cl: double;
  ran: TVec3d;
begin

  randseed:= 1;

  AssignFile(myFile, 'Segments.txt');
  ReWrite(myFile);

  Writeln(myFile, '$MeshFormat');
  Writeln(myFile, '2.1 0 8');
  Writeln(myFile, '$EndMeshFormat');
  Writeln(myFile, '$Nodes');
  Writeln(myFile, 2*length(segments)+length(points));

  for k:= 0 to length(segments)-1 do begin

    kk:= GET_INDEX(segments[k].name, false);
  end;
  for k:= 0 to length(points)-1 do begin

    kk:= GET_INDEX(points[k].name, true);
  end;
  for k:= 0 to length(segments)-1 do begin

    ran:= GET_RAND(segments[k].name, false);
    if segments[k].namei <> -1 then begin

      ran:= GET_RAND(segments[k].namei, true);
    end;
    Writeln(myFile, 2*k, ' ',
            segments[k].nodei[0]+ran[0], ' ',
            segments[k].nodei[1]+ran[1], ' ',
            segments[k].nodei[2]+ran[2]
            );
    if segments[k].namej <> -1 then begin

      ran:= GET_RAND(segments[k].namej, true);
    end;
    Writeln(myFile, 2*k+1, ' ',
            segments[k].nodej[0]+ran[0], ' ',
            segments[k].nodej[1]+ran[1], ' ',
            segments[k].nodej[2]+ran[2]
            );
  end;
  for k:= 0 to length(points)-1 do begin

    ran:= GET_RAND(points[k].name, true);
    Writeln(myFile, 2*length(segments)+k, ' ',
            points[k].nodei[0]+ran[0], ' ',
            points[k].nodei[1]+ran[1], ' ',
            points[k].nodei[2]+ran[2]
            );
  end;
  Writeln(myFile, '$EndNodes');
  Writeln(myFile, '$Elements');
  Writeln(myFile, length(segments)+length(points));
  for k:= 0 to length(segments)-1 do begin

    Writeln(myFile, k, ' 1 2 0 ',
            segments[k].name, ' ',
            2*k, ' ',
            2*k+1
            );
  end;
  for k:= 0 to length(points)-1 do begin

    Writeln(myFile, length(segments)+k, ' 15 2 0 ',
            points[k].name, ' ',
            2*length(segments)+k
            );
  end;
  Writeln(myFile, '$EndElements');
  Writeln(myFile, '$ElementData');
  Writeln(myFile, '1');
  Writeln(myFile, '"explicit"');
  Writeln(myFile, '1');
  Writeln(myFile, '0.0');
  Writeln(myFile, '3');
  Writeln(myFile, '0');
  Writeln(myFile, '1');
  if length(state) > 0 then begin

    kk:= state[0].first+state[0].firstP;
  end else begin

    kk:= length(segments)+length(points);
  end;
  Writeln(myFile, kk);
  for k:= 0 to kk-1 do begin

    Writeln(myFile, k, ' 0.0');
  end;
  Writeln(myFile,'$EndElementData');

  for kk:= 0 to length(state)-1 do begin

    Writeln(myFile, '$ElementData');
    Writeln(myFile, '1');
    Writeln(myFile, '"implicit"');
    Writeln(myFile, '1');
    Writeln(myFile, '0.0');
    Writeln(myFile, '3');
    Writeln(myFile, kk);
    Writeln(myFile, '1');
    Writeln(myFile, state[kk].last-state[kk].first+1+state[kk].lastP-state[kk].firstP+1);
    for k:= state[kk].first to state[kk].last do begin

      cl:= GET_COLOR(segments[k].name, false);
      Writeln(myFile, k, ' ', cl);
    end;
    for k:= state[kk].firstP to state[kk].lastP do begin

      cl:= GET_COLOR(points[k].name, true);
      Writeln(myFile, length(segments)+k, ' ', cl);
    end;
    Writeln(myFile,'$EndElementData');
  end;
  CloseFile(myFile);
end;
//==============================================================================

function POINT_LINE_PROJECTION(const point: TVec3d;
                               const r0: TVec3d;
                               const rt: TVec3d
                               ): double;
begin

  result:= DOT(DIF(point, r0), rt)/DOT(rt, rt);
end;
//==============================================================================

function POINT_CURVE_PROJECTION(const point: TVec3d;
                                const curvePoint: TParametricCurvePoint;
                                const param: TParam;
                                const tol: double = 1e-10
                                ): double;
var

  P: TParametricCurvePoint;
  b: boolean;
begin

  P:= curvePoint;
  b:= true;
  result:= P.t;
  while b do begin

    result:= result+POINT_LINE_PROJECTION(point, P.r, P.rt);
    b:= abs(result-P.t)*NORM(P.rt) > tol;
    result:= param.CONVERT(result);
    P:= param.EVALUATE_FUNCTION_AND_DERIVATIVES(result);
  end;
end;
//==============================================================================

procedure INTERSECT_LINE_LINE(const r0_t, r0_s: TVec3d;
                              const rt, rs: TVec3d;
                              var t, s: double);
var

  rt_norm, rs_norm: double;
  a, b, c, d, value, dota: double;
  dx: TVec3d;
begin

  rt_norm:= sqrt(DOT(rt, rt));
  rs_norm:= sqrt(DOT(rs, rs));
  dota:= DOT(rt, rs);
  dx:= DIF(r0_s, r0_t);
  //Solving the 2x2 system of linear eq
  a:= rt_norm*rt_norm;
  b:= -dota;
  c:= dota;
  d:= -rs_norm*rs_norm;
  value:= 1/(a*d-b*c);

  //s and t of aparent intersection
  t:= value*(d*DOT(dx, rt) - b*DOT(dx, rs));
  s:= value*(-c*DOT(dx, rt) + a*DOT(dx, rs));
end;
//
//dota -> DOT(rt, rs)
//dx -> r0_s - r0_t
//rt -> vel of "t" (diff)
//==============================================================================

function PROJECT_BISECTOR_LINES(
  const r0_t, r0_s: TVec3d;
  const rt, rs: TVec3d
  ): double;
var

  direction, n: TVec3d;
  value: double;
begin

  n:= CROSS(rt, rs);
  value:= DOT(n, n);
  // check if lines are non-parallel
  if value > 1e-100 then begin

    value:= SQRT(DOT(rs, rs)/DOT(rt, rt));
    if DOT(rt, rs) < 0 then begin

      value:= -value;
    end;
    // bisector
    direction[0]:= rt[0]*value + rs[0];
    direction[1]:= rt[1]*value + rs[1];
    direction[2]:= rt[2]*value + rs[2];
    // perpendicular to bisector
    direction:= CROSS(direction, n);
    INTERSECT_LINE_LINE(r0_t, r0_s, direction, rs, value, result);
  end else begin

    result:= DOT(DIF(r0_t, r0_s), rs)/DOT(rs, rs);
  end;
end;
//
// returns the s coordinate of the projected r0_t
//==============================================================================

function PROJECT_BISECTOR(
  const point1, point2: TParametricCurvePoint;
  const Eval2: ParamFunc1D;
  const WithinTol2: FuncAreWithinTol;
  const tol: double = 1e-10
  ): TParametricCurvePoint;
var

  t0, value: double;
begin

  result:= point2;
  repeat begin

    t0:= result.t;
    value:= PROJECT_BISECTOR_LINES(point1.r, result.r, point1.rt, result.rt);
    result:= Eval2(t0+value);
  end;
  until WithinTol2(result.t, t0, tol);
end;
//
// returns the projection of point1 on curve 2
//==============================================================================

function FIND_CURVES_DEBONDING(
  const Eval1, Eval2: ParamFunc1D;
  const point1, point2: TParametricCurvePoint;
  const tol: double;
  const inward: boolean
  ): TParametricCurvePoint;
var

  k: integer;
  p1, p2: TParametricCurvePoint;
  a, b, c, d, value, D1, D2: double;
  n, aa, dx: TVec3d;
begin

  p1:= point1;
  p2:= point2;

  k:= 0;
  while k <= 5 do begin

    CROSS(p1.rt, p2.rt, n);
    CROSS(n, p2.rt, aa);
    dx:= DIF(p2.r, p1.r);
    if inward then begin

      value:= tol/NORM(aa);
    end else begin

      value:= -tol/NORM(aa);
    end;
    dx[0]:= dx[0]+value*aa[0];
    dx[1]:= dx[1]+value*aa[1];
    dx[2]:= dx[2]+value*aa[2];
    a:= DOT(p1.rt, p2.rt);
    b:= -DOT(p2.rt, p2.rt);
    c:= DOT(p1.rt, p1.rt);
    d:= -DOT(p1.rt, p2.rt);
    value:=a*d-b*c;
    if abs(value) > 1e-100 then begin

      value:= 1/value;
      D1:= value*(d*DOT(dx, p2.rt)-b*DOT(dx, p1.rt));
      D2:= value*(-c*DOT(dx, p2.rt)+a*DOT(dx, p1.rt));
      p1:= Eval1(p1.t+D1);
      p2:= Eval2(p2.t+D2);
    end else begin

      p1.t:= 1e-308;
      if inward then begin

        p1.t:= -p1.t;
      end;
      k:= 1000000;
    end;
    k:= k+1;
  end;

  result:= p1;
end;
//
//  - find parameter t of the first curve which distance to the second curve is tol
//  - distance is measured perpendicular to curve 2 (direction of aa)
//  - dx is the vector from the current parapoint in curve 1 to the current parapoint in curve 2
//
//==============================================================================

function IS_INSIDE(
  const t: double;
  const t1, t2: double;
  const inverted: boolean = false
  ): boolean;
var

  ti, tj: double;
begin

  if not inverted then begin

    ti:= t1;
    tj:= t2;
  end else begin

    ti:= t2;
    tj:= t1;
  end;
  if tj > ti then begin

    // normal case
    result:= (t >= ti) and (t <= tj);
  end else begin

    // in case of parameterization jump (e.g. circle)
    result:= (t >= ti) or (t <= tj);
  end;
end;
//
//  check if parameter t is strickly inside (ignoring tolerances) interval [ti, tj]
//  ti, tj - parameters of the start (inward = true) and end (inward = false) points
//  inverted - true if ti and tj are not in the natural order (acording to previous line)
//
//==============================================================================

function OVERLAP(
  const bi1, bj1: TParametricCurvePoint;                 // boundaries of obj 1
  const bi2, bj2: TParametricCurvePoint;                 // boundaries of obj 1
  const debondi1, debondj1: TParametricCurvePoint;       // debonding points of obj 1
  const debondi2, debondj2: TParametricCurvePoint;       // debonding points of obj 2
  const inverted: boolean;                               // parameterizations point in opposite directions
  const tol: double;                                     // tol for point projection
  const Eval1, Eval2: ParamFunc1D;
  var min1, max1, min2, max2: double;                    // overlap initial and final points on both curves
  const WithinTol1: FuncAreWithinTol;
  const WithinTol2: FuncAreWithinTol
  ): boolean;                                            // true if overlap exists, false otherwise
var

  value: double;
begin

  // check if an overlap between debonding and boundary intervals exist (otherwise no overlap possible)
  result:= false;
  if IS_INSIDE(bi2.t, debondi2.t, debondj2.t, inverted) or IS_INSIDE(bj2.t, debondi2.t, debondj2.t, inverted) or
     IS_INSIDE(debondi2.t, bi2.t, bj2.t) then begin

    if IS_INSIDE(bi1.t, debondi1.t, debondj1.t) or IS_INSIDE(bj1.t, debondi1.t, debondj1.t) or
       IS_INSIDE(debondi1.t, bi1.t, bj1.t) then begin

      // min is the ovelap start param of curve 1
      if IS_INSIDE(bi1.t, debondi1.t, debondj1.t) then begin

        min1:= bi1.t;
        min2:= PROJECT_BISECTOR(bi1, debondi2, Eval2, WithinTol2, 0.001*tol).t;
      end else begin

        min1:= debondi1.t;
        min2:= debondi2.t;
      end;
      if not inverted then begin

        if IS_INSIDE(bi2.t, min2, debondj2.t) then begin

          min1:= PROJECT_BISECTOR(bi2, debondi1, Eval1, WithinTol1, 0.001*tol).t;
          min2:= bi2.t;
        end;
      end else begin

        if IS_INSIDE(bj2.t, debondj2.t, min2) then begin

          min1:= PROJECT_BISECTOR(bj2, debondi1, Eval1, WithinTol1, 0.001*tol).t;
          min2:= bj2.t;
        end;
      end;

      // max is the ovelap end param of curve 1
      if IS_INSIDE(bj1.t, debondi1.t, debondj1.t) then begin

        max1:= bj1.t;
        max2:= PROJECT_BISECTOR(bj1, debondi2, Eval2, WithinTol2, 0.001*tol).t;
      end else begin

        max1:= debondj1.t;
        max2:= debondj2.t;
      end;
      if not inverted then begin

        if IS_INSIDE(bj2.t, debondi2.t, max2) then begin

          max1:= PROJECT_BISECTOR(bj2, debondi1, Eval1, WithinTol1, 0.001*tol).t;
          max2:= bj2.t;
        end;
      end else begin

        if IS_INSIDE(bi2.t, max2, debondi2.t) then begin

          max1:= PROJECT_BISECTOR(bi2, debondi1, Eval1, WithinTol1, 0.001*tol).t;
          max2:= bi2.t;
        end;
      end;

      // check if min and max are within a tol
      if WithinTol1(min1, max1, tol) or IS_INSIDE(max1, min1, debondj1.t) then begin

        result:= true;
      end;
    end;
  end;
end;
//==============================================================================

function PointProjections(const point: TVec3d;
                          const Curve: TStraightLine;
                          const direction: TVec3d):double;
var
    dx: TVec3d;
    rt_NormSq: double;
    t, s: double;
    r_average: TVec3d;
    rt_norm, rs_norm: double;
    dota: double;

begin

  rt_norm:= sqrt(DOT(Curve.rt0, Curve.rt0));
  rs_norm:= sqrt(DOT(direction, direction));
  dota:= DOT(Curve.rt0, direction);

  dx:= DIF(point, Curve.r0);
  INTERSECT_LINE_LINE(curve.r0, point, curve.rt0, direction, t, s);

  Result:= t;

end;
//Projects a point using a direction vector
//==============================================================================

function ParallelDirection(const Param_t: TStraightLine): TVec3d;

begin

  if abs(Param_t.rt0[0]) > abs(Param_t.rt0[1]) then begin

    if abs(Param_t.rt0[0]) > abs(Param_t.rt0[2]) then begin

      Result[0]:= Param_t.rt0[1];
      Result[1]:= -Param_t.rt0[0];
      Result[2]:= 0;

    end else begin

      Result[0]:= -Param_t.rt0[2];
      Result[1]:= 0;
      Result[2]:= Param_t.rt0[0];

    end;

  end else begin

    if abs(Param_t.rt0[1]) > abs(Param_t.rt0[2]) then begin

      Result[0]:= -Param_t.rt0[1];
      Result[1]:= Param_t.rt0[0];
      Result[2]:= 0;

    end else begin

      Result[0]:= -Param_t.rt0[2];
      Result[1]:= 0;
      Result[2]:= Param_t.rt0[0];

    end;

  end;


end;
//THIS IS PERPENDICULAR DIRECTION!!!!
//==============================================================================

function INTERSECT_SEGMENT_SEGMENT(
  const obj1, obj2: TObject3D;
  const tol: double;
  const tolL: double;
  const tolA: double
  ): TObjectArray;
var

  dist, s, t, ti, tj, si, sj, dd, tolH, value, rt_norm, rs_norm: double;
  n, dx: TVec3d;
  Param_t, Param_s: TStraightLine;
  p: pointer;
  sine: double;

  dota: double;

  min, max: double;
  proj1, proj2: double;

  debug: double;

  direction: TVec3d; //from t to s
  r_average: TVec3d;
  a: double;

  k: integer;
  x: array[0..3] of double;


begin

    //Assign pointers to the parametric information
    p:= obj1.Param;
    Param_t:= p;
    p:= obj2.Param;
    Param_s:= p;

    CROSS(Param_t.rt0, Param_s.rt0, n);
    dota:= DOT(Param_t.rt0, Param_s.rt0);
    rt_norm:= sqrt(DOT(Param_t.rt0, Param_t.rt0));
    rs_norm:= sqrt(DOT(Param_s.rt0, Param_s.rt0));
    value:= sqrt(DOT(n, n));
    sine:= value/(rt_norm*rs_norm);

    if dota > 0 then begin

      r_average[0]:= 0.5*(Param_t.rt0[0]/rt_norm + Param_s.rt0[0]/rs_norm);
      r_average[1]:= 0.5*(Param_t.rt0[1]/rt_norm + Param_s.rt0[1]/rs_norm);
      r_average[2]:= 0.5*(Param_t.rt0[2]/rt_norm + Param_s.rt0[2]/rs_norm);

    end else begin

      r_average[0]:= 0.5*(Param_t.rt0[0]/rt_norm - Param_s.rt0[0]/rs_norm);
      r_average[1]:= 0.5*(Param_t.rt0[1]/rt_norm - Param_s.rt0[1]/rs_norm);
      r_average[2]:= 0.5*(Param_t.rt0[2]/rt_norm - Param_s.rt0[2]/rs_norm);

    end;

    //Get boundaries of curve 1
    if obj1.Inward[0] = True then begin

      ti:= obj1.ParaPoint[0];
      tj:= obj1.ParaPoint[1];

    end else begin

      tj:= obj1.ParaPoint[0];
      ti:= obj1.ParaPoint[1];

    end;

    //Get boundaries of curve 2
    if obj2.Inward[0] = True then begin

      si:= obj2.ParaPoint[0];
      sj:= obj2.ParaPoint[1];

    end else begin

      sj:= obj2.ParaPoint[0];
      si:= obj2.ParaPoint[1];

    end;

    if value >= 1e-10 then begin //Not parallel

      // Normalizing the n vector
      n[0]:= n[0]/value;
      n[1]:= n[1]/value;
      n[2]:= n[2]/value;

      dx:= DIF(Param_s.r0, Param_t.r0);
      dist:= DOT(dx, n); //Dist between parallel planes

      CROSS(r_average, n, direction);

      if abs(dist) <= tol then begin

        INTERSECT_LINE_LINE(Param_t.r0, Param_s.r0, Param_t.rt0, Param_s.rt0, t, s);

        //distance between the aparent intersection and the debounding points
        dd:= sqrt(tol*tol - dist*dist) / sine;

        if (rt_norm*ABS(t-ti) <= dd) or (rt_norm*ABS(t-tj) <= dd) or
           ((ti-t)*(tj-t) < 0) then begin

          if (rs_norm*ABS(s-si) <= dd) or (rs_norm*ABS(s-sj) <= dd) or
             ((si-s)*(sj-s) < 0) then begin

             //Calculate min
             min:= t-dd/rt_norm;
             if ti > min then begin min:= ti; end;
             if dota > 0 then begin

               value:= si;

             end else begin

               value:= sj;
             end;
             value:= PointProjections(Param_s.r(value), Param_t, direction);
             if value > min then begin

              min:= value;

             end;


             //Calculate max
             max:= t+dd/rt_norm;
             if tj < max then begin max:= tj; end;
             if dota > 0 then begin

               value:= sj;

             end else begin

               value:= si;
             end;
             value:= PointProjections(Param_s.r(value), Param_t, direction);
             if value < max then begin

              max:= value;

             end;


             if (min-max)*rt_norm <= tol then begin

                //Assign New Results
                SetLength(Result, 3);
                Result[0]:= TPoint.Create;
                Result[1]:= TPoint.Create;
                Result[2]:= TCurve.Create;
                SetLength(Result[2].Boundary, 2);
                Result[2].Boundary[0]:= Result[0];
                Result[2].Boundary[1]:= Result[1];

                //To avoid Curves with min > max. Creates a Curve with 0 length on both ends
                x[0]:= min;
                x[1]:= max;
                x[2]:= min;
                x[3]:= max;
                if min > max then begin
                  if min = ti then begin

                    x[1]:= min;
                    x[2]:= max;
                  end else if max = tj then begin

                    x[0]:= max;
                    x[3]:= min;
                  end;
                end;

                //Assign newContains
                k:= Length(obj1.Contain);
                SetLength(obj1.Contain,k+1);
                obj1.Contain[k].Obj:= Result[2];
                SetLength(obj1.Contain[k].ParaPoint,2);
                SetLength(obj1.Contain[k].Inward,2);
                obj1.Contain[k].ParaPoint[0]:= x[0];
                obj1.Contain[k].ParaPoint[1]:= x[1];
                obj1.Contain[k].Inward[0]:= True;
                obj1.Contain[k].Inward[1]:= False;

                //Assign newContains (OBJ2)
                k:= Length(obj2.Contain);
                SetLength(obj2.Contain,k+1);
                obj2.Contain[k].Obj:= Result[2];
                SetLength(obj2.Contain[k].ParaPoint,2);
                SetLength(obj2.Contain[k].Inward,2);
                obj2.Contain[k].ParaPoint[0]:= PointProjections(Param_t.r(x[2]), Param_s, direction);
                obj2.Contain[k].ParaPoint[1]:= PointProjections(Param_t.r(x[3]), Param_s, direction);

                //Inwards change order if dot product between 2 rt0 is -1
                if dota > 0 then begin

                  obj2.Contain[k].Inward[0]:= True;
                  obj2.Contain[k].Inward[1]:= False;

                end else begin

                  obj2.Contain[k].Inward[0]:= False;
                  obj2.Contain[k].Inward[1]:= True;

                end;

             end;

          end;

        end;

      end;

    end else begin //Parallel

       direction:= ParallelDirection(Param_t);

       value:= PointProjections(Param_t.r(obj1.ParaPoint[0]), Param_s, direction);

       dd:= sqrt(POWER(Param_t.r(obj1.ParaPoint[0])[0] - Param_s.r(value)[0], 2) +
            POWER(Param_t.r(obj1.ParaPoint[0])[1] - Param_s.r(value)[1], 2) +
            POWER(Param_t.r(obj1.ParaPoint[0])[2] - Param_s.r(value)[2], 2));

       if dd < tol then
       begin


             min:= ti;
             if dota > 0 then begin

               value:= si;

             end else begin

               value:= sj;
             end;

             proj1:= PointProjections(Param_s.r(value), Param_t, direction);
             if proj1 > min then begin min:= proj1; end;

             max:= tj;
             if dota > 0 then begin

               value:= sj;

             end else begin

               value:= si;
             end;

             proj2:= PointProjections(Param_s.r(value), Param_t, direction);
             if proj2 < max then begin max:= proj2; end;

             if (min-max)*rt_norm <= tol then begin

                //Assign New Results
                SetLength(Result, 3);
                Result[0]:= TPoint.Create;
                Result[1]:= TPoint.Create;
                Result[2]:= TCurve.Create;
                SetLength(Result[2].Boundary, 2);
                Result[2].Boundary[0]:= Result[0];
                Result[2].Boundary[1]:= Result[1];

                //To avoid Curves with min > max. Creates a Curve with 0 length on both ends
                x[0]:= min;
                x[1]:= max;
                x[2]:= min;
                x[3]:= max;
                if min > max then begin
                  if min = ti then begin

                    x[1]:= min;
                    x[2]:= max;
                  end else if max = tj then begin

                    x[0]:= max;
                    x[3]:= min;
                  end;
                end;

                //Assign newContains
                k:= Length(obj1.Contain);
                SetLength(obj1.Contain,k+1);
                obj1.Contain[k].Obj:= Result[2];
                SetLength(obj1.Contain[k].ParaPoint,2);
                SetLength(obj1.Contain[k].Inward,2);
                obj1.Contain[k].ParaPoint[0]:= x[0];
                obj1.Contain[k].ParaPoint[1]:= x[1];
                obj1.Contain[k].Inward[0]:= True;
                obj1.Contain[k].Inward[1]:= False;

                //Assign newContains
                k:= Length(obj2.Contain);
                SetLength(obj2.Contain,k+1);
                obj2.Contain[k].Obj:= Result[2];
                SetLength(obj2.Contain[k].ParaPoint,2);
                SetLength(obj2.Contain[k].Inward,2);
                obj2.Contain[k].ParaPoint[0]:= PointProjections(Param_t.r(x[2]), Param_s, direction);
                obj2.Contain[k].ParaPoint[1]:= PointProjections(Param_t.r(x[3]), Param_s, direction);

                //Inwards change order if dot product between 2 rt0 is -1
                if dota > 0 then begin

                  obj2.Contain[k].Inward[0]:= True;
                  obj2.Contain[k].Inward[1]:= False;

                end else begin

                  obj2.Contain[k].Inward[0]:= False;
                  obj2.Contain[k].Inward[1]:= True;

                end;

             end;

       end;

    end;

end;
//==============================================================================

function INTERSECT_POINT_POINT(
  const obj1, obj2: TObject3D;
  const tol: double = 1e-3;
  const tolL: double = 1e-3;
  const tolA: double = 1e-3
  ): TObjectArray;
var

  k: integer;
  Param1, Param2: TParam0D;
  p: pointer;
begin

    p:= obj1.Param;
    Param1:= p;
    p:= obj2.Param;
    Param2:= p;
    if sqrt(POWER(Param2.r[0] - Param1.r[0], 2) +
            POWER(Param2.r[1] - Param1.r[1], 2) +
            POWER(Param2.r[2] - Param1.r[2], 2)) < tol then
    begin

      setLength(Result, 1);
      Result[0]:= TPoint.Create;
      k:= Length(obj1.Contain);
      SetLength(obj1.Contain, k+1);
      obj1.Contain[k].Obj:= Result[0];
      k:= Length(obj2.Contain);
      SetLength(obj2.Contain, k+1);
      obj2.Contain[k].Obj:= Result[0];
    end else begin

      setLength(Result, 0);
    end;
end;
//==============================================================================

function INTERSECT_POINT_SEGMENT(
  const obj1, obj2: TObject3D;
  const tol: double = 1e-3;
  const tolL: double = 1e-3;
  const tolA: double = 1e-3
  ): TObjectArray;
var

  k: integer;
  Param1: TParam0D;
  Param3: TParam1D;
  p: pointer;
  dx: TVec3d;
  rt_NormSq, ti, tj, t: double;
begin
    //Assign pointers to the parametric information
    p:= obj1.Param;
    Param1:= p;
    p:= obj2.Param;
    Param3:= p;
    dx:= DIF(Param1.r, Param3.r(0));
    rt_NormSq:= DOT(Param3.rt(0), Param3.rt(0));
    t:= DOT(dx, Param3.rt(0))/rt_NormSq;
    if DOT(dx,dx)-rt_NormSq*t*t <= tol*tol then begin

      if obj2.Inward[0] = True then begin

        ti:= obj2.ParaPoint[0];
        tj:= obj2.ParaPoint[1];
      end else begin

        tj:= obj2.ParaPoint[0];
        ti:= obj2.ParaPoint[1];
      end;
      if t < ti then begin

        if intpower(ti-t,2)*rt_NormSq <= intpower(tol,2) then begin

          setLength(Result, 1);
        end;
      end else begin

        if t > tj then begin

          if intpower(t-tj,2)*rt_NormSq <= intpower(tol,2) then begin

            setLength(Result, 1);
          end;
        end else begin

          setLength(Result, 1);
        end;
      end;
      if Length(Result) = 1 then begin

        Result[0]:= TPoint.Create;
        Result[0].Param:= nil;

        k:= Length(obj1.Contain);
        SetLength(obj1.Contain, k+1);
        obj1.Contain[k].Obj:= Result[0];
        k:= Length(obj2.Contain);
        SetLength(obj2.Contain, k+1);
        obj2.Contain[k].Obj:= Result[0];
        SetLength(obj2.Contain[k].ParaPoint, 1);
        obj2.Contain[k].ParaPoint[0]:= t;
      end;
    end;
end;
//==============================================================================

function INTERSECT_POINT_CURVE(
  const obj1, obj2: TObject3D;
  const tol: double = 1e-3;
  const tolL: double = 1e-3;
  const tolA: double = 1e-3
  ): TObjectArray;
var

  kkk: integer;
  Pair: TClusterPairs;
  Curve: TCurve;
  p: pointer;
  point: TParametricCurvePoint;
  t: double;

begin

    obj1.GENERATE_BASE_CLUSTERS(0);

    // safe step from parameterization
    obj2.GENERATE_BASE_CLUSTERS(obj2.Param.SAFE_STEP);

    FIND_PAIRS(0, length(obj2.Cluster)-1,
               obj1.Cluster, obj2.Cluster,
               tol,
               Pair);

    if length(Pair) <> 0 then begin

      p:= obj2;
      Curve:= p;

      point:= obj2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES(Curve.net[Pair[0].index2].p[0]);       //  Curve.net[Pair[0].index2].point
      t:= POINT_CURVE_PROJECTION(obj1.Param.EVALUATE, point, Curve.Param, 1e-3*tol);

      // distance from point to its projection
      if NORM(DIF(obj1.Param.EVALUATE, obj2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES(t).r)) <= tol then begin

        if Curve.IS_INSIDE(t, tol) then begin

          setlength(result, 1);
          result[0]:= TPoint.Create;
          kkk:= length(obj1.Contain);
          setlength(obj1.Contain, kkk+1);
          obj1.Contain[kkk].Obj:= result[0];
          kkk:= length(obj2.Contain);
          setlength(obj2.Contain, kkk+1);
          obj2.Contain[kkk].Obj:= result[0];
          setlength(obj2.Contain[kkk].ParaPoint, 1);
          obj2.Contain[kkk].ParaPoint[0]:= t;
        end;
      end;
    end;
end;
//==============================================================================

function INTERSECT_CURVE_CURVE(
  const obj1, obj2: TObject3D;
  const tol: double = 1e-3;
  const tolL: double = 1e-3;
  const tolA: double = 1e-3
  ): TObjectArray;
var

  k, kk, kkk: integer;
  p: pointer;
  Curve1, Curve2: TCurve;
  bol: boolean;
  groups: TGroupsOfPairs;
  low, high: integer;

  debondi1, debondj1, debondi2, debondj2: TParametricCurvePoint;
  bi1, bj1, bi2, bj2: TParametricCurvePoint;
  min1, max1, min2, max2: double;
begin


    if obj1.Cluster[0].radius < obj2.Cluster[0].radius then begin

      p:= obj1;
      Curve1:= p;
      p:= obj2;
      Curve2:= p;
    end else begin

      p:= obj2;
      Curve1:= p;
      p:= obj1;
      Curve2:= p;
    end;


    groups:= FIND_GROUPS(Curve1.Cluster, Curve2.Cluster,
                         Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                         Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                         tol,
                         Curve1.net, Curve2.net,
                         Curve1.Param.AVERAGE_T,
                         Curve1.Param.AVERAGE_T,
                         Curve2.Param.AVERAGE_T,
                         Curve2.Param.AVERAGE_T);
    // sort neighbors of curve 1 net nodes in increasing order
    kk:= 0;
    for k:= 1 to length(Curve1.net)-1 do begin

      kkk:= Curve1.net[kk].neighbor[1];
      if Curve1.net[kkk].neighbor[0] <> kk then begin

         Curve1.net[kkk].neighbor[1]:= Curve1.net[kkk].neighbor[0];
         Curve1.net[kkk].neighbor[0]:= kk;
      end;
      kk:= kkk;
    end;


    // handle each group
    for k:= 0 to groups.NumGroups-1 do begin


      // find pair index with lowest (inward = true) and highest (inward = false) parameter t
      for kkk:= 0 to length(Curve1.net)-1 do begin

        if (Curve1.net[kkk].order = k) then begin

          if (Curve1.net[kkk].neighbor[0] = -1) or
             ((Curve1.net[kkk].neighbor[0] <> -1) and
              (Curve1.net[Curve1.net[kkk].neighbor[0]].order <> k)) then begin

            low:= Curve1.net[kkk].pair[0];
          end;
          if (Curve1.net[kkk].neighbor[1] = -1) or
             ((Curve1.net[kkk].neighbor[1] <> -1) and
              (Curve1.net[Curve1.net[kkk].neighbor[1]].order <> k)) then begin

            high:= Curve1.net[kkk].pair[0];
          end;
        end;
      end;

      debondi1:= FIND_CURVES_DEBONDING(
                   Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve1.net[groups.Pair[low].index1].point,
                   Curve2.net[groups.Pair[low].index2].point,
                   tol, true);
      debondj1:= FIND_CURVES_DEBONDING(
                   Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve1.net[groups.Pair[high].index1].point,
                   Curve2.net[groups.Pair[high].index2].point,
                   tol, false);
      bol:= DOT(Curve1.net[groups.Pair[low].index1].rt,
                Curve2.net[groups.Pair[low].index2].rt) > 0;
      debondi2:= FIND_CURVES_DEBONDING(
                   Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve2.net[groups.Pair[low].index2].point,
                   Curve1.net[groups.Pair[low].index1].point,
                   tol, bol);
      debondj2:= FIND_CURVES_DEBONDING(
                   Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                   Curve2.net[groups.Pair[high].index2].point,
                   Curve1.net[groups.Pair[high].index1].point,
                   tol, not bol);

      if not Curve1.isLoop then begin

        bi1:= Curve1.startPoint;
        bj1:= Curve1.endPoint;
      end else begin

        bi1:= debondi1;
        bj1:= debondj1;
      end;

      if not Curve2.isLoop then begin

        bi2:= Curve2.startPoint;
        bj2:= Curve2.endPoint;
      end else begin

        bi2:= debondi2;
        bj2:= debondj2;
      end;


      if OVERLAP(bi1, bj1,
                 bi2, bj2,
                 debondi1, debondj1,
                 debondi2, debondj2,
                 not bol, tol,
                 Curve1.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                 Curve2.Param.EVALUATE_FUNCTION_AND_DERIVATIVES,
                 min1, max1,
                 min2, max2,
                 Curve1.Param.ARE_WITHIN_TOL,
                 Curve2.Param.ARE_WITHIN_TOL
                 ) then begin

        // create objects
        kk:= length(result);
        setlength(result, kk+3);
        result[kk]:= TPoint.Create;
        result[kk+1]:= TPoint.Create;
        result[kk+2]:= TCurve.Create;
        setlength(result[kk+2].Boundary, 2);
        result[kk+2].Boundary[0]:= result[kk];
        result[kk+2].Boundary[1]:= result[kk+1];

        // obj1
        kkk:= length(Curve1.Contain);
        setlength(Curve1.Contain, kkk+1);
        Curve1.Contain[kkk].Obj:= result[kk+2];
        setlength(Curve1.Contain[kkk].ParaPoint, 2);
        Curve1.Contain[kkk].ParaPoint[0]:= min1;
        Curve1.Contain[kkk].ParaPoint[1]:= max1;
        setlength(Curve1.Contain[kkk].Inward, 2);
        Curve1.Contain[kkk].Inward[0]:= true;
        Curve1.Contain[kkk].Inward[1]:= false;

        // obj2
        kkk:= length(Curve2.Contain);
        setlength(Curve2.Contain, kkk+1);
        Curve2.Contain[kkk].Obj:= result[kk+2];
        setlength(Curve2.Contain[kkk].ParaPoint, 2);
        Curve2.Contain[kkk].ParaPoint[0]:= min2;
        Curve2.Contain[kkk].ParaPoint[1]:= max2;
        setlength(Curve2.Contain[kkk].Inward, 2);
        Curve2.Contain[kkk].Inward[0]:= bol;
        Curve2.Contain[kkk].Inward[1]:= not bol;
      end;
    end;

end;
//==============================================================================

function INTERSECT(const obj1, obj2: TObject3D;
                   const tol: double;
                   const tolL: double;
                   const tolA: double
                   ): TObjectArray;
var

  k: integer;
begin

  setLength(Result, 0);


  // Point - Point intersection
  if (obj1.ObjType = otPoint) and
     (obj2.ObjType = otPoint) then begin

    Result:= INTERSECT_POINT_POINT(obj1, obj2, tol, tolL, tolA);
  end;


  // Point - curve intersection
  if (obj1.ObjType = otPoint) and
     (obj2.ObjType = otCurve) then begin

    if obj2.Param.ParamType = ptStraightLine then begin


      Result:= INTERSECT_POINT_SEGMENT(obj1, obj2, tol, tolL, tolA);
    end else begin

      Result:= INTERSECT_POINT_CURVE(obj1, obj2, tol, tolL, tolA);
    end;
  end;


  {// Point - Surface intersection
  if ((obj1.ObjType = otPoint) and (obj2.ObjType = otSurface)) or
     ((obj1.ObjType = otSurface) and (obj2.ObjType = otPoint)) then begin
    if obj1.ObjType = otPoint then begin

      p:=obj1;
      Point1:= p;
      p:=obj2;
      Surf1:= p;
    end else begin

      p:=obj2;
      Point1:= p;
      p:=obj1;
      Surf1:= p;
    end;

    // Point - Plane intersection
    if Surf1.Param.ParamType = ptPlane then begin


    end;

  end;        }




  // Curve - Curve
  if (obj1.ObjType = otCurve) and
     (obj2.ObjType = otCurve) then begin

    if (obj1.Param.ParamType = ptStraightLine) and
       (obj2.Param.ParamType = ptStraightLine) then begin


      Result:= INTERSECT_SEGMENT_SEGMENT(obj1, obj2, tol, tolL, tolA);
    end else begin


      Result:= INTERSECT_CURVE_CURVE(obj1, obj2, tol, tolL, tolA);
    end;
  end;

  // add contained by to all result objects
  for k:=0 to Length(Result)-1 do begin

    SetLength(Result[k].ContainedBy, 2);
    Result[k].ContainedBy[0]:= obj1;
    Result[k].ContainedBy[1]:= obj2;
  end;
end;
//==============================================================================
 {
function TNetNode1D.t: double;
begin

  result:= self.point.t;
end;
//==============================================================================

function TNetNode1D.r0: TVec3d;
begin

  result:= self.point.r0;
end;
//==============================================================================

function TNetNode1D.rt: TVec3d;
begin

  result:= self.point.rt;
end;
//==============================================================================     }
function RECONCILE(const obj: TObjectArray;
                     const tol:double=1e-3;
                     const tolL:double=1e-3;
                     const tolA:double=1e-3
                     ): TObjectArray; Overload;

  //============================================================================

  procedure QuickSort(var List: TProtoObjectArray; iLo, iHi: integer) ;
  var
    Lo       : integer;
    Hi       : integer;
    T_obj    : TObject3d;
    T_para   : Double;
    T_inwards: boolean;
    Mid      : Double;

  begin
    Lo := iLo;
    Hi := iHi;
    Mid:= List[(Lo + Hi) div 2].ParaPoint[0];
    repeat

      while List[Lo].ParaPoint[0] < Mid do Inc(Lo) ;
      while List[Hi].ParaPoint[0] > Mid do Dec(Hi) ;

      if Lo <= Hi then
      begin
        T_obj                 := List[Lo].Obj;
        T_para                := List[Lo].ParaPoint[0];
        T_inwards             := List[Lo].Inward[0];

        List[Lo].ParaPoint[0] := List[Hi].ParaPoint[0];
        List[Lo].Obj          := List[Hi].Obj;
        List[Lo].Inward[0]    := List[Hi].Inward[0];
        List[Hi].ParaPoint[0] := T_para;
        List[Hi].Obj          := T_obj;
        List[Hi].Inward[0]    := T_inwards;
        Inc(Lo);
        Dec(Hi);
      end;

    until Lo > Hi;

    if Hi > iLo then QuickSort(List, iLo, Hi);
    if Lo < iHi then QuickSort(List, Lo, iHi);

  end;
  //============================================================================

  function FIND_CONTAINED_INDEX(const obj1, curve: TObject3d): integer;
  var
    i: integer;

  begin

    Result:= -1;
    for i:=0 to Length(obj1.Contain)-1 do begin

      if obj1.Contain[i].Obj = curve then begin

        Result:= i;
      end;
    end;
  end;
  //============================================================================

  function ContainedCurves(const Curve, actualObject: TObject3d):TObjectArray;
  var
    i, numContained: integer;

  begin

    numContained:= Length(Curve.ContainedBy);
    SetLength(Result, 0);
    for i:=0 to numContained-1 do begin

      if actualObject <> Curve.ContainedBy[i] then begin
        SetLength(Result, Length(Result)+1);
        Result[Length(Result)-1]:= Curve.ContainedBy[i];
      end;
    end;
  end;
  //Generare Array with all the objects that are contained by an active curve
  //============================================================================

  procedure BreakPoint(const t: double;
                       const Curve: TObject3d;
                       const newCurve: TObject3d;
                       const Obj: TObject3d;
                       const trueIndex: integer);

  var

    k, j: integer;
  begin

    k:= length(Obj.Contain);
    setlength(Obj.Contain, k+1);
    Obj.Contain[k].obj:= newCurve;
    setlength(Obj.Contain[k].ParaPoint, 2);
    setlength(Obj.Contain[k].Inward, 2);

    j:= FIND_CONTAINED_INDEX(Obj, Curve);

    Obj.Contain[k].ParaPoint[0]:= Obj.Contain[j].ParaPoint[trueIndex];
    Obj.Contain[k].ParaPoint[1]:= t;                                  //Change t to projection
    Obj.Contain[j].ParaPoint[trueIndex]:= t;

    if trueIndex = 0 then begin

      Obj.Contain[k].Inward[0]:= Obj.Contain[j].Inward[0];
      Obj.Contain[k].Inward[1]:= Obj.Contain[j].Inward[1];
    end else begin

      Obj.Contain[k].Inward[0]:= Obj.Contain[j].Inward[1];
      Obj.Contain[k].Inward[1]:= Obj.Contain[j].Inward[0];
    end;
  end;
  //Change the object (break it into smaller pieces)
  //1. Add new curve between the break point
  //2. Change the original explicit and implicit objects so they have this point and curve
  //============================================================================

  procedure SORT_EQUALS(var ProtoArray: TProtoObjectArray);

  var
    i, j, k, kk, kkk, numObjs: integer;
    current: TProtoObject3d;
    numEquals: integer;
    ProtoObj: TProtoObject3d;
    bol: boolean;

  begin

    numObjs:= Length(ProtoArray);

    for i:=0 to numObjs-2 do begin

      if Abs(ProtoArray[i].ParaPoint[0] - ProtoArray[i+1].ParaPoint[0]) < 0.0001 then begin

        //Find NumEquals (begins at "i" and ends at "k"
        k:= i+2;
        while (k < Length(ProtoArray)) and
              (Abs(ProtoArray[k].ParaPoint[0] - ProtoArray[k-1].ParaPoint[0]) < 0.0001) do begin

          inc(k, 1);
        end;

        inc(k, -1);

        bol:= False;
        kk:=0;
        while (kk < i) do begin

          for kkk:=i to k do begin

            if ProtoArray[kkk].Obj = ProtoArray[kk].Obj then begin

              ProtoObj:= ProtoArray[kkk];
              ProtoArray[kkk]:= ProtoArray[i];
              ProtoArray[i]:= ProtoObj;

              kk:=i; //Breaks the kk while

              bol:= True;
            end;
          end;

          inc(kk, 1);
        end;

        if not bol then begin

          kk:= i;
          while (kk <= k) and (ProtoArray[kk].Inward[0] = True) do begin

            inc(kk, 1);
          end;

          if kk <> k+1 then begin

            for kkk:=i to k do begin

              if (ProtoArray[kkk].Obj = ProtoArray[kk].Obj) and (ProtoArray[kkk].Inward[0] = True) then begin

                 ProtoObj:= ProtoArray[kkk];
                 ProtoArray[kkk]:= ProtoArray[i];
                 ProtoArray[i]:= ProtoObj;
              end;
            end;

          end else begin

            kk:= k+1;
            while (kk < Length(ProtoArray)) do begin

              for kkk:=i to k do begin

                if ProtoArray[kkk].Obj = ProtoArray[kk].Obj then begin

                  ProtoObj:= ProtoArray[kkk];
                  ProtoArray[kkk]:= ProtoArray[i];
                  ProtoArray[i]:= ProtoObj;

                  kk:=Length(ProtoArray); //Breaks the kk while
                end;
              end;

              inc(kk, 1);
            end;
          end;
        end;
      end;
    end;
  end;
  //Swap objs, if t is the same
  //============================================================================

  procedure ADD_TO_LIST(const obj: TObject3d; var ObjArray: TObjectArray);
  var

    i: integer;
    bol: boolean;
  begin

    bol:= false;
    for i:= 0 to length(ObjArray)-1 do begin

      if objArray[i] = Obj then begin

        bol:= true;
      end;
    end;

    if bol = false then begin

      setlength(ObjArray, length(ObjArray)+1);
      ObjArray[length(ObjArray)-1]:= obj;
    end;
  end;
  //Adds a new object to the Results
  //============================================================================

  procedure ELIMINATE_ACTIVE(const obj: TObject3d; var Active: TObjectArray);
  var

    i, j: integer;
  begin

    i:= 0;
    while (Active[i] <> obj) do begin

      inc(i);
    end;

    for j:= i+1 to length(Active)-1 do begin

      Active[j-1]:= Active[j];
    end;

    SetLength(Active, length(Active)-1);
  end;
  //============================================================================

  procedure ELIMINATE_NILS(const obj: TObject3d);
  var

    j, k: integer;
  begin

    j:= 0;
    while j < length(obj.Contain) do begin

      if obj.Contain[j].Obj = nil then begin

        for k:= j to length(obj.Contain)-2 do begin

          obj.Contain[k]:= obj.Contain[k+1];
        end;

        setlength(obj.Contain, length(obj.Contain)-1);
      end else begin

      j:= j+1;
      end;
    end;
  end;
  //============================================================================

  function CHECK_CHAIN(obj1, obj2: TProtoObject3d): boolean;
  begin

    result:= false;
    if (obj1.Obj.Boundary[1] = obj2.Obj.Boundary[0]) then begin

       result:= true;
    end;
  end;
  //============================================================================

  function INWARD(const contain: TProtoObject3d; const bol: boolean): integer;
  begin

    if contain.Inward[0] = bol then begin

      result:= 0;
    end else begin

      result:= 1;
    end;
  end;
  //Returns the index of the parapoint (0 or 1) of some INWARD
  //============================================================================

  function FIND_OTHER_CONTAINEDBY(const objs: TObjectArray; const obj: TObject3d): TObject3d;
  begin

    if obj = objs[0] then begin

      Result:= objs[1];
    end else begin

      Result:= objs[0];
    end;
  end;
  //============================================================================

  function FIND_PARAMETER_POINT(
    const point, obj: TObject3d
    ): double;
  var

    i: integer;
    value, min, max: double;
  begin

    min:= Infinity;
    max:= -Infinity;
    for i:= 0 to length(obj.Contain)-1 do begin

      if obj.Contain[i].obj.ObjType = otCurve then begin

        if obj.Contain[i].obj.boundary[0] = point then begin

          value:= obj.Contain[i].Parapoint[0];
          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
        if obj.Contain[i].obj.boundary[1] = point then begin

          value:= obj.Contain[i].Parapoint[1];
          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
      end;
    end;
    if (min = Infinity) or (min = max) then begin

      result:= min;
    end else begin

      result:= 0.5*(min+max);
    end;
  end;
  //============================================================================

  function FIND_PARAMETER_POINT2(
    const point, obj: TObject3d;
    const groups: TGroups
    ): double;
  var

    i: integer;
    value, min, max: double;
  begin

    min:= Infinity;
    max:= -Infinity;
    for i:= 0 to length(obj.Contain)-1 do begin

      if obj.Contain[i].obj.ObjType = otCurve then begin

        if (obj.Contain[i].obj.boundary[0] = point) or
           groups.ARE_ASSOCIATED(obj.Contain[i].obj.boundary[0], point) then begin

          value:= obj.Contain[i].Parapoint[0];
          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
        if (obj.Contain[i].obj.boundary[1] = point) or
           groups.ARE_ASSOCIATED(obj.Contain[i].obj.boundary[1], point) then begin

          value:= obj.Contain[i].Parapoint[1];
          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
      end;
    end;
    if (min = Infinity) or (min = max) then begin

      result:= min;
    end else begin

      result:= 0.5*(min+max);
    end;
  end;
  //============================================================================
    {
  function FIND_PROJECTED_PARAMETER(
    const point, obj1, obj3: TObject3d;
    var points: TProtoObjectArray
    ): double;
  var

    i, j, k, count: integer;
    t1, t2: double;
    x: TVec3d;
    obj2: TObject3d;
    bol: boolean;
    p: pointer;
    Param0D: TParam0D;
    Param1D: TParam1D;
  begin

    count:= 0;
    result:= 0;
    for j:= 0 to length(obj1.Contain)-1 do begin

      if obj1.Contain[j].Obj.ObjType = otPoint then begin

        if obj1.Contain[j].Obj = point then begin

          p:= obj1.Param;
          Param1D:= p;
          t1:= obj1.Contain[j].ParaPoint[0];
          x:= Param1D.EVALUATE_FUNCTION_AND_DERIVATIVES(t1);
          p:= obj3.Param;
          Param1D:= p;
          result:= result+POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0));
          inc(count);
        end;
      end else begin

        for k:= 0 to 1 do begin

          if obj1.Contain[j].Obj.Boundary[k] = point then begin

            bol:= false;
            for i:= 0 to length(points)-1 do begin

              if point = points[i].Obj then begin

                x[0]:= points[i].ParaPoint[0];
                x[1]:= points[i].ParaPoint[1];
                x[2]:= points[i].ParaPoint[2];
                bol:= true;
              end;
            end;
            if bol = false then begin

              obj2:= FIND_OTHER_CONTAINEDBY(obj1.Contain[j].obj.ContainedBy, obj1);
              t2:= FIND_PARAMETER_POINT(point, obj2);
              t1:= obj1.Contain[j].ParaPoint[k];
              p:= obj2.Param;
              Param1D:= p;
              x2:= Param1D.r(t2);
              p:= obj1.Param;
              Param1D:= p;
              x1:= Param1D.r(t1);
              rt1:= Param1D.rt(0);
              x[0]:= 0.5*(x1[0]+x2[0]);
              x[1]:= 0.5*(x1[1]+x2[1]);
              x[2]:= 0.5*(x1[2]+x2[2]);

              kk:= length(points);
              setlength(points, kk+1);
              points[kk].Obj:= point;
              setlength(points[kk].ParaPoint, 3);
              points[kk].ParaPoint[0]:= x[0];
              points[kk].ParaPoint[1]:= x[1];
              points[kk].ParaPoint[2]:= x[2];
            end;

            //Project
            p:= obj3.Param;
            Param1D:= p;
            result:= POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0)); //Projection might not be possible of x

          end;
        end;
      end;
    end;
  end;                 }
  //
  //  points contains the average x, y, z coordinates of the point dual definition
  //
  //============================================================================
      {
  procedure ADD_POINT_COORDS(
    const point: TObject3d;
    const curve: TObject3d;
    var points: TProtoObjectArray);
  var

    i, j, k, kk: integer;
    t1, t2: double;
    x: TVec3d;

    obj1, obj2: TObject3d;

    bol: boolean;
    p: pointer;
    Param0D: TParam0D;
    Param1D: TParam1D;
  begin

    // check if point is already in the list of points
    bol:= false;
    for i:= 0 to length(points)-1 do begin

      if point = points[i].Obj then begin

        bol:= true;
      end;
    end;

    if not bol then begin

      if curve = nil then begin

        // a point directly contained in a obj

      end else begin

        if point = then


      end;
    end;



    for i:= 0 to length(obj)-1 do begin

      for j:= 0 to length(obj[i].Contain)-1 do begin

        if obj[i].Contain[j].Obj.ObjType = otPoint then begin

          if obj[i].Contain[j].Obj = point then begin

            p:= obj[i].Param;
            Param1D:= p;
            t1:= obj[i].Contain[j].ParaPoint[0];
            x:= Param1D.EVALUATE_FUNCTION_AND_DERIVATIVES(t1);
            p:= obj3.Param;
            Param1D:= p;
            result:= result+POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0));
            inc(count);
          end;
        end else begin

          for k:= 0 to 1 do begin

            if obj[i].Contain[j].Obj.Boundary[k] = point then begin

              if bol = false then begin

                obj2:= FIND_OTHER_CONTAINEDBY(obj1.Contain[j].obj.ContainedBy, obj1);
                t2:= FIND_PARAMETER_POINT(point, obj2);
                t1:= obj[i].Contain[j].ParaPoint[k];
                p:= obj2.Param;
                Param1D:= p;
                x2:= Param1D.r(t2);
                p:= obj[i].Param;
                Param1D:= p;
                x1:= Param1D.r(t1);
                rt1:= Param1D.rt(0);
                x[0]:= 0.5*(x1[0]+x2[0]);
                x[1]:= 0.5*(x1[1]+x2[1]);
                x[2]:= 0.5*(x1[2]+x2[2]);

                kk:= length(points);
                setlength(points, kk+1);
                points[kk].Obj:= point;
                setlength(points[kk].ParaPoint, 3);
                points[kk].ParaPoint[0]:= x[0];
                points[kk].ParaPoint[1]:= x[1];
                points[kk].ParaPoint[2]:= x[2];
              end;

              //Project
              p:= obj3.Param;
              Param1D:= p;
              result:= POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0)); //Projection might not be possible of x

            end;
          end;
        end;
      end;
    end;
  end;       }
  //
  //  points contains the average x, y, z coordinates of the point dual definition
  //
  //============================================================================
    {
  function FIND_POINT(
    const point: TObject3d;
    const x: TVec3d;
    var list: TProtoObjectArray
    ): TVec3d;
  var
    i, kk: integer;
    bol: boolean;

  begin

     bol:= false;
     for i:= 0 to length(list)-1 do begin

       if point = list[i].Obj then begin

         result[0]:= list[i].ParaPoint[0];
         result[1]:= list[i].ParaPoint[1];
         result[2]:= list[i].ParaPoint[2];

         bol:= true;
       end;
     end;

     if bol = false then begin

       kk:= length(list);
       setlength(list, kk+1);
       list[kk].Obj:= point;
       setlength(list[kk].ParaPoint, 3);
       list[kk].ParaPoint[0]:= x[0];
       list[kk].ParaPoint[1]:= x[1];
       list[kk].ParaPoint[2]:= x[2];

       result:= x;
     end;
  end;              }
  //============================================================================

  function FIND_BREAK_PARAM(
    const obj1, obj3: TObject3d;
    const Break: TProtoObject3d;
    const isInitiation: boolean;
    const groups: TGroups;
    var points: TProtoObjectArray
    ): double;
  var

    i, kk, kkk, index2: integer;
    bol: boolean;
    p: pointer;
    Param1D: TParam1D;
    x2, x1 ,rt1, x: TVec3d;
    obj2, point: TObject3d;
    t1, t2: double;
    value, min, max: double;
  begin

    index2:= FIND_CONTAINED_INDEX(obj1, Break.Obj);
    kkk:= INWARD(obj1.Contain[index2], isInitiation);
    point:= Break.obj.boundary[kkk];
    kkk:= groups.GROUP_INDEX(point);
    bol:= false;
    if kkk = -1 then begin

      result:= FIND_PARAMETER_POINT(point, obj3);
      if True then

    end else begin
       {
      for i:= 0 to length(groups.group[kkk])-1 do begin

        result:= FIND_PARAMETER_POINT(groups.group[kkk][i], obj3);
      end;
       }
      bol:= false;
      min:= Infinity;
      max:= -Infinity;
      for i:= 0 to length(groups.group[kkk])-1 do begin

        value:= FIND_PARAMETER_POINT(groups.group[kkk][i], obj3);
        if value <> Infinity then begin

          bol:= true;
          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
      end;
      if (bol = false) or (min =  max) then begin

        result:= min;
      end else begin

        result:= 0.5*(min+max);
      end;

    end;


    if not bol then begin

      bol:= false;
      for i:= 0 to length(points)-1 do begin

        if point = points[i].Obj then begin

          x[0]:= points[i].ParaPoint[0];
          x[1]:= points[i].ParaPoint[1];
          x[2]:= points[i].ParaPoint[2];
          bol:= true;
        end;
      end;

      if bol = false then begin

        obj2:= FIND_OTHER_CONTAINEDBY(Break.obj.ContainedBy, obj1);
        t2:= FIND_PARAMETER_POINT(point, obj2);
        t1:= Break.ParaPoint[0];
        p:= obj2.Param;
        Param1D:= p;
        x2:= Param1D.r(t2);
        p:= obj1.Param;
        Param1D:= p;
        x1:= Param1D.r(t1);
        rt1:= Param1D.rt(0);
        x[0]:= 0.5*(x1[0]+x2[0]);
        x[1]:= 0.5*(x1[1]+x2[1]);
        x[2]:= 0.5*(x1[2]+x2[2]);

        kk:= length(points);
        setlength(points, kk+1);
        points[kk].Obj:= point;
        setlength(points[kk].ParaPoint, 3);
        points[kk].ParaPoint[0]:= x[0];
        points[kk].ParaPoint[1]:= x[1];
        points[kk].ParaPoint[2]:= x[2];
      end;

      //Project
      p:= obj3.Param;
      Param1D:= p;
      result:= POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0)); //Projection might not be possible of x
    end;
  end;
  //============================================================================
   {
  function FIND_BREAK_PARAM2(
    const obj1, obj3: TObject3d;
    const point: TObject3d;
    const groups: TGroups;
    var points: TProtoObjectArray
    ): double;
  var

    i, j, kk, kkk, index2: integer;
    bol: boolean;
    p: pointer;
    Param1D: TParam1D;
    x2, x1 ,rt1, x: TVec3d;
    obj2: TObject3d;
    t1, t2: double;

    value, min, max: double;
  begin

    kkk:= groups.GROUP_INDEX(point);
    if kkk = -1 then begin

      result:= FIND_PARAMETER_POINT(point, obj3);    // unlikely that result is different from infinite
    end else begin

      min:= Infinite;
      max:= -Infinite;
      for i:= 0 to length(groups.group[kkk])-1 do begin

        value:= FIND_PARAMETER_POINT(groups.group[kkk][i], obj3);
        if not IsInfinite(value) then begin

          if value < min then begin

            min:= value;
          end;
          if value > max then begin

            max:= value;
          end;
        end;
      end;
      if IsInfinite(min) or (min =  max) then begin

        result:= min;
      end else begin

        result:= 0.5*(min+max);
      end;
    end;


    if IsInfinite(result) then begin

      for j:= 0 to length(obj1.Contain)-1 do begin

        if obj1.Contain[j].Obj.ObjType = otCurve then begin

          if groups.GROUP_INDEX(obj1.Contain[j].Obj.Boundary[0]) = kkk then begin

            min;
          end;
          if groups.GROUP_INDEX(obj1.Contain[j].Obj.Boundary[1]) = kkk then begin

            min;
          end;

        end;
      end;


      bol:= false;
      for i:= 0 to length(points)-1 do begin

        if point = points[i].Obj then begin

          x[0]:= points[i].ParaPoint[0];
          x[1]:= points[i].ParaPoint[1];
          x[2]:= points[i].ParaPoint[2];
          bol:= true;
        end;
      end;

      if bol = false then begin

        obj2:= FIND_OTHER_CONTAINEDBY(Break.obj.ContainedBy, obj1);
        t2:= FIND_PARAMETER_POINT(point, obj2);
        t1:= Break.ParaPoint[0];
        p:= obj2.Param;
        Param1D:= p;
        x2:= Param1D.r(t2);
        p:= obj1.Param;
        Param1D:= p;
        x1:= Param1D.r(t1);
        rt1:= Param1D.rt(0);
        x[0]:= 0.5*(x1[0]+x2[0]);
        x[1]:= 0.5*(x1[1]+x2[1]);
        x[2]:= 0.5*(x1[2]+x2[2]);

        kk:= length(points);
        setlength(points, kk+1);
        points[kk].Obj:= point;
        setlength(points[kk].ParaPoint, 3);
        points[kk].ParaPoint[0]:= x[0];
        points[kk].ParaPoint[1]:= x[1];
        points[kk].ParaPoint[2]:= x[2];
      end;

      //Project
      p:= obj3.Param;
      Param1D:= p;
      result:= POINT_LINE_PROJECTION(x, Param1D.r(0), Param1D.rt(0));
    end;
  end;                  }
  //============================================================================

  function CONTAINED_BY_ARE_EQUAL(
    const list1: TObjectArray;
    const list2: TObjectArray
    ): boolean;
  var

    k, kk, kkk: integer;
    bol: boolean;
  begin

    kkk:= length(list1);
    if length(list2) = kkk then begin

      result:= true;
      for k:= 0 to kkk-1 do begin

        bol:= false;
        for kk:= 0 to kkk-1 do begin

          if list1[k] = list2[kk] then begin

            bol:= true;
          end;
        end;
        if not bol then begin

          result:= false;
        end;
      end;
    end else begin

      result:= false;
    end;
  end;
  //============================================================================
   {
  function PAIR(
    const obj: TObject3d;
    const pairs: TPairs
    ): TObject3d;
  var

    k: integer;
  begin

    result:= nil;
    for k:= 0 to length(pairs)-1 do begin

      if pairs[k][0] = obj then begin

        result:= pairs[k][1];
      end;
      if pairs[k][1] = obj then begin

        result:= pairs[k][0];
      end;
    end;
  end;          }
  //============================================================================
var

  i, j, k, kk, kkk, kkkk, numObjs, numProtoObjs, numContains, numActives: integer;
  Params: array of TStraightLine;
  ProtoObjects: TProtoObjectArray;
  points: TProtoObjectArray;
  currentActive: TObject3d;
  currentProtoObj: TProtoObject3d;
  max, min, ti, tb, tf: double;
  Active, ActiveContains, other: TObjectArray;
  ActiveCopy: TObject3d;
  newCurve: TObject3d;
  index1, index2, index3:integer;

  obj1, obj2, obj3: TObject3d;
  done: boolean;
  bol, bol1, bol2: boolean;
  p: pointer;
  Param1D: TParam1D;
  Chain: TObjectArray;
  value: double;
  toCollapse: TObjectArray;
  a,b: integer;
  point: TObject3d;

  r0, rt, x1, x2: TVec3d;
  groups: TGroups;
  gmsh: TGmshFile;
  g1, g2, g3, g4: integer;
  t1, t2, t3, t4: double;
  created: TObjectArray;
  rt0: double;
  b0, b1: TObject3d;

begin

  gmsh.ADD_EXPLICIT(obj);
  gmsh.ADD_STATE(obj);


  // create list of created objects
  numObjs:= length(obj);
  setlength(created, 0);
  for i:= 0 to numObjs-1 do begin

    for j:= 0 to length(obj[i].Contain)-1 do begin

      ADD_TO_LIST(obj[i].Contain[j].Obj, created);
      if obj[i].Contain[j].Obj.ObjType = otCurve then begin

        ADD_TO_LIST(obj[i].Contain[j].Obj.Boundary[0], created);
        ADD_TO_LIST(obj[i].Contain[j].Obj.Boundary[1], created);
      end;
    end;
  end;

  // group points that are within tolL
  groups.ASSOCIATION(obj, tolL);

// reduce
//groups.REDUCE(obj);


  done:= false;
  setlength(points, 0);
  while not done do begin

    done:= true;
    for i:= 0 to numObjs-1 do begin

      numContains:= length(obj[i].Contain);
      if (numContains > 0) and (obj[i].ObjType = otCurve) then begin

        //Get the parameters of obj[i] to use in BreakPoint
        p:= obj[i].Param;
        Param1D:= p;

        r0:= Param1D.r(0);
        rt:= Param1D.rt(0);

        numProtoObjs:= 0;
        SetLength(ProtoObjects, 0);

        for j:= 0 to numContains-1 do begin

          //If curve
          if obj[i].Contain[j].Obj.ObjType = otCurve then begin

            numProtoObjs:= numProtoObjs + 2;
            SetLength(ProtoObjects, numProtoObjs);

            ProtoObjects[numProtoObjs-2].Obj:= obj[i].Contain[j].Obj;
            ProtoObjects[numProtoObjs-1].Obj:= obj[i].Contain[j].Obj;

            SetLength(ProtoObjects[numProtoObjs-2].ParaPoint, 1);
            SetLength(ProtoObjects[numProtoObjs-1].ParaPoint, 1);
            ProtoObjects[numProtoObjs-2].ParaPoint[0]:= obj[i].Contain[j].ParaPoint[0];
            ProtoObjects[numProtoObjs-1].ParaPoint[0]:= obj[i].Contain[j].ParaPoint[1];

            SetLength(ProtoObjects[numProtoObjs-2].Inward, 1);
            SetLength(ProtoObjects[numProtoObjs-1].Inward, 1);
            ProtoObjects[numProtoObjs-2].Inward[0]:= obj[i].Contain[j].Inward[0];
            ProtoObjects[numProtoObjs-1].Inward[0]:= obj[i].Contain[j].Inward[1];

          end;
        end;

        //Order this list and do breakPoints each time there are overlapping objs
        QuickSort(ProtoObjects, 0, Length(ProtoObjects)-1);
        SORT_EQUALS(ProtoObjects);
        numProtoObjs:= length(ProtoObjects);
        setlength(Active, 0);
        for j:= 0 to numProtoObjs-1 do begin

          // INITIATION
          if (ProtoObjects[j].Inward[0] = true) and (length(Active) > 0) then begin

            for k:= 0 to length(Active)-1 do begin

              index1:= FIND_CONTAINED_INDEX(Obj[i], Active[k]);
              index2:= FIND_CONTAINED_INDEX(Obj[i], ProtoObjects[j].Obj);
              if (ProtoObjects[j].ParaPoint[0] <> obj[i].Contain[index1].Parapoint[0]) and
                 (ProtoObjects[j].ParaPoint[0] <> obj[i].Contain[index1].Parapoint[1]) and
                 (not groups.ARE_ASSOCIATED(ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(true)],
                                            Active[k].Boundary[0])) and
                 (not groups.ARE_ASSOCIATED(ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(true)],
                                            Active[k].Boundary[1])) then begin



                done:= false;
                obj3:= FIND_OTHER_CONTAINEDBY(obj[i].Contain[index1].obj.ContainedBy, obj[i]);



                tb:= FIND_BREAK_PARAM(obj[i], obj3, ProtoObjects[j], true, groups, points);

                point:= ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(true)];
                //t1:= FIND_PARAMETER_POINT2(point, obj[i], groups);

                kkk:= obj[i].Contain[index1].INWARD_INDEX(true);

                //Create New Curve
                newCurve:= TCurve.Create;
                ADD_TO_LIST(newCurve, created);
                setlength(newCurve.Boundary, 2);
                newCurve.Boundary[0]:= Active[k].Boundary[kkk];                 // may be different boundary!
                newCurve.Boundary[1]:= point;
                setlength(newCurve.ContainedBy, 2);
                newCurve.ContainedBy[0]:= Active[k].ContainedBy[0];
                newCurve.ContainedBy[1]:= Active[k].ContainedBy[1];
                Active[k].Boundary[kkk]:= point;


                  for kk:= 0 to 1 do begin

                    if Active[k].ContainedBy[kk] = obj[i] then begin

                      value:= ProtoObjects[j].ParaPoint[0];
                    end else begin

                      value:= tb;
                    end;

                    BreakPoint(value,
                            Active[k],
                            newCurve,
                            Active[k].ContainedBy[kk],
                            kkk);
                  end;


                gmsh.ADD_STATE(obj);
              end;
            end;
            setlength(Active, length(Active)+1);
            Active[length(Active)-1]:= ProtoObjects[j].Obj;

          end;

          // FINALIZATION
          if (ProtoObjects[j].Inward[0] = false) and (length(Active) > 1) then begin

            ELIMINATE_ACTIVE(ProtoObjects[j].obj, Active);

            for k:= 0 to length(Active)-1 do begin

              index1:= FIND_CONTAINED_INDEX(Obj[i], Active[k]);
              index2:= FIND_CONTAINED_INDEX(Obj[i], ProtoObjects[j].Obj);
              if (ProtoObjects[j].ParaPoint[0] <> obj[i].Contain[index1].Parapoint[1]) and
                 (ProtoObjects[j].ParaPoint[0] <> obj[i].Contain[index1].Parapoint[0]) and
                 (not groups.ARE_ASSOCIATED(ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(false)],
                                           Active[k].Boundary[0])) and
                 (not groups.ARE_ASSOCIATED(ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(false)],
                                           Active[k].Boundary[1])) then begin

                done:= false;
                obj3:= FIND_OTHER_CONTAINEDBY(obj[i].Contain[index1].obj.ContainedBy, obj[i]);

                tb:= FIND_BREAK_PARAM(obj[i], obj3, ProtoObjects[j], false, groups, points);

                point:= ProtoObjects[j].Obj.Boundary[obj[i].Contain[index2].INWARD_INDEX(false)];
              //  t1:= FIND_PARAMETER_POINT2(point, obj[i], groups);

                kkk:= obj[i].Contain[index1].INWARD_INDEX(true);

                //Create New Curve
                newCurve:= TCurve.Create;
                ADD_TO_LIST(newCurve, created);
                setlength(newCurve.Boundary, 2);
                newCurve.Boundary[0]:= Active[k].Boundary[kkk];
                newCurve.Boundary[1]:= point;
                setlength(newCurve.ContainedBy, 2);
                newCurve.ContainedBy[0]:= Active[k].ContainedBy[0];
                newCurve.ContainedBy[1]:= Active[k].ContainedBy[1];
                Active[k].Boundary[kkk]:= point;


                for kk:= 0 to 1 do begin

                  if Active[k].ContainedBy[kk] = obj[i] then begin

                    value:= ProtoObjects[j].ParaPoint[0];
                  end else begin

                    value:= tb;
                  end;

                  BreakPoint(value,
                           Active[k],
                           newCurve,
                           Active[k].ContainedBy[kk],
                           kkk);
                end;

                gmsh.ADD_STATE(obj);
              end;
            end;
          end;

          // Also for the first active
          if length(Active) = 0 then begin

            setlength(Active, 1);
            Active[0]:= ProtoObjects[j].Obj;

          end else if (ProtoObjects[j].Inward[0] = False) and
                      (length(Active) = 1) and
                      (ProtoObjects[j].obj = Active[0]) then begin

            setlength(Active, 0);
          end;
        end;
      end;
    end;
  end;


  // Union
  for i:= 0 to numObjs-1 do begin

    numContains:= length(obj[i].Contain);
    if (numContains > 0) and (obj[i].ObjType = otCurve) then begin

      for j:= 0 to numContains-1 do begin

        if (obj[i].Contain[j].Obj <> nil) and (obj[i].Contain[j].Obj.ObjType = otCurve) then begin

          for k:= j+1 to numContains-1 do begin

            if (obj[i].Contain[k].Obj <> nil) and (obj[i].Contain[k].Obj.ObjType = otCurve) then begin

              g1:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[0]);
              if g1 = -1 then begin

                g1:= Integer(Pointer(obj[i].Contain[j].Obj.Boundary[0]))
              end;
              g2:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[1]);
              if g2 = -1 then begin

                g2:= Integer(Pointer(obj[i].Contain[j].Obj.Boundary[1]))
              end;
              g3:= groups.GROUP_INDEX(obj[i].Contain[k].Obj.Boundary[0]);
              if g3 = -1 then begin

                g3:= Integer(Pointer(obj[i].Contain[k].Obj.Boundary[0]))
              end;
              g4:= groups.GROUP_INDEX(obj[i].Contain[k].Obj.Boundary[1]);
              if g4 = -1 then begin

                g4:= Integer(Pointer(obj[i].Contain[k].Obj.Boundary[1]))
              end;
              // compare boundaries
              if ((g1 = g3) or
                  (g1 = g4)) and
                 ((g2 = g3) or
                  (g2 = g4)) and
                  (g1 <> g2) then begin

                kkk:= length(obj[i].Contain[k].Obj.ContainedBy);
                obj3:= obj[i].Contain[k].Obj;

                for kk:= 0 to kkk-1 do begin

                  obj2:= obj3.ContainedBy[kk];
                  obj[i].Contain[j].Obj.ADD_CONTAINED_BY(obj2);
                  index1:= FIND_CONTAINED_INDEX(obj2, obj[i].Contain[j].Obj);
                  index2:= FIND_CONTAINED_INDEX(obj2, obj3);
                  if index1 = -1 then begin

                    obj2.Contain[index2].Obj:= obj[i].Contain[j].Obj;
                  end else begin

                    g1:= obj2.Contain[index1].INWARD_INDEX(true);
                    t1:= obj2.Contain[index1].ParaPoint[g1];
                    g2:= obj2.Contain[index1].INWARD_INDEX(false);
                    t2:= obj2.Contain[index1].ParaPoint[g2];
                    g3:= obj2.Contain[index2].INWARD_INDEX(true);
                    t3:= obj2.Contain[index2].ParaPoint[g1];
                    g4:= obj2.Contain[index2].INWARD_INDEX(false);
                    t4:= obj2.Contain[index2].ParaPoint[g2];
                    if t3 < t1 then begin

                      obj2.Contain[index1].ParaPoint[g1]:= t3;
                    end;
                    if t4 > t2 then begin

                      obj2.Contain[index1].ParaPoint[g2]:= t4;
                    end;
                    obj2.Contain[index2].Obj:= nil;
                  end;

                  gmsh.ADD_STATE(obj);
                end;
              end;
            end;
          end;
        end;
      end;

      //Resize "j" list
      ELIMINATE_NILS(obj[i]);
    end;
  end;



  // new collapse
  for k:= 0 to length(groups.group)-1 do begin

    for i:= 0 to numObjs-1 do begin

      if (obj[i].ObjType = otCurve) and (obj[i].Param <> nil) then begin

        rt0:= NORM(obj[i].Param.EVALUATE_FUNCTION_AND_DERIVATIVES(0).rt);
      end;

      bol:= false;
      min:= Infinity;
      max:= -Infinity;
      for j:= 0 to length(obj[i].Contain)-1 do begin

        if obj[i].Contain[j].Obj.ObjType = otPoint then begin

          kk:= groups.GROUP_INDEX(obj[i].Contain[j].Obj);
          if kk = k then begin

            obj[i].Contain[j].Obj:= nil;
            if obj[i].Contain[j].ParaPoint[0] < min then begin

              min:= obj[i].Contain[j].ParaPoint[0];
            end;
            if obj[i].Contain[j].ParaPoint[0] > max then begin

              max:= obj[i].Contain[j].ParaPoint[0];
            end;
          end;
        end else begin

          kk:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[0]);
          kkk:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[1]);
          if kk = k then begin

            if obj[i].Contain[j].ParaPoint[0] < min then begin

              min:= obj[i].Contain[j].ParaPoint[0];
            end;
            if obj[i].Contain[j].ParaPoint[0] > max then begin

              max:= obj[i].Contain[j].ParaPoint[0];
            end;
          end;
          if kkk = k then begin

            if obj[i].Contain[j].ParaPoint[1] < min then begin

              min:= obj[i].Contain[j].ParaPoint[1];
            end;
            if obj[i].Contain[j].ParaPoint[1] > max then begin

              max:= obj[i].Contain[j].ParaPoint[1];
            end;
          end;
          if (kk = k) and (kkk = k) then begin

            // collapse curve
            obj[i].Contain[j].Obj:= nil;
          end;
          if (kk = k) xor (kkk = k) then begin

            bol:= true;  // point is contained in a non-collapsing curve
            if kk = k then begin

              obj[i].Contain[j].Obj.Boundary[0]:= groups.group[k][0];
            end else begin

              obj[i].Contain[j].Obj.Boundary[1]:= groups.group[k][0];
            end;
          end;
        end;
      end;

      if min <> Infinity then begin

        bol1:= false;
        bol2:= false;
        if  abs(max-obj[i].ParaPoint[obj[i].INWARD_INDEX(false)])*rt0 < tolL then begin

          value:= obj[i].ParaPoint[obj[i].INWARD_INDEX(false)];
          bol1:= true;
        end;
        if abs(min-obj[i].ParaPoint[obj[i].INWARD_INDEX(true)])*rt0 < tolL then begin

          value:= obj[i].ParaPoint[obj[i].INWARD_INDEX(true)];
          bol2:= true;
        end;

        if bol1 = bol2 then begin

          value:= 0.5*(min+max);
        end;

        if not bol then begin

          kk:= length(obj[i].Contain);
          setlength(obj[i].Contain, kk+1);
          obj[i].Contain[kk].Obj:= groups.group[k][0];
          setlength(obj[i].Contain[kk].ParaPoint, 1);
          obj[i].Contain[kk].ParaPoint[0]:= value;
        end else begin

          for j:= 0 to length(obj[i].Contain)-1 do begin

            if (obj[i].Contain[j].Obj <> nil) and (obj[i].Contain[j].Obj.ObjType = otCurve)  then begin

              kk:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[0]);
              kkk:= groups.GROUP_INDEX(obj[i].Contain[j].Obj.Boundary[1]);
              if k = kk then begin

                obj[i].Contain[j].ParaPoint[0]:= value;
              end;
              if k = kkk then begin

                obj[i].Contain[j].ParaPoint[1]:= value;
              end;
            end;
          end;
        end;
      end;
      ELIMINATE_NILS(Obj[i]);
    end;
  end;



  // add implicit objects to result
  for i:= 0 to numObjs-1 do begin

    for j:= 0 to length(obj[i].Contain)-1 do begin

      if obj[i].Contain[j].Obj.ObjType = otCurve then begin

        ADD_TO_LIST(obj[i].Contain[j].Obj.Boundary[0], result);
        ADD_TO_LIST(obj[i].Contain[j].Obj.Boundary[1], result);
        kk:= length(obj[i].Contain[j].Obj.ContainedBy);
        for k:= 0 to kk-1 do begin

          obj[i].Contain[j].Obj.Boundary[0].ADD_CONTAINED_BY(obj[i].Contain[j].Obj.ContainedBy[k]);
          obj[i].Contain[j].Obj.Boundary[1].ADD_CONTAINED_BY(obj[i].Contain[j].Obj.ContainedBy[k]);
        end;
      end else begin

        ADD_TO_LIST(obj[i].Contain[j].Obj, result);
      end;
    end;
  end;
  for i:= 0 to numObjs-1 do begin

    for j:= 0 to length(obj[i].Contain)-1 do begin

      ADD_TO_LIST(obj[i].Contain[j].Obj, result);
    end;
  end;


  // free lost objects
  for i:= 0 to length(created)-1 do begin

    bol:= false;
    for j:= 0 to length(result)-1 do begin

      if created[i] = result[j] then begin

        bol:= true;
      end;
    end;
    if not bol then begin

      created[i].Free;
    end;
  end;


  gmsh.ADD_STATE(obj);
  gmsh.SAVE_TO_FILE;
end;
//==============================================================================

constructor TPoint.Create;
begin

  self.ObjType:= otPoint;
end;
//==============================================================================

constructor TCurve.Create;
begin

  self.ObjType:= otCurve;
end;
//==============================================================================

constructor TSurface.Create;
begin

  self.ObjType:= otSurface;
end;
//==============================================================================

constructor TSolid.Create;
begin

  self.ObjType:= otSolid;
end;
//==============================================================================
//==============================================================================

procedure TPoint.GENERATE_BASE_CLUSTERS(const safeStep: double);
begin

  setlength(Cluster, 1);
  Cluster[0].clusterType:= ctSphere;
  Cluster[0].center:= self.Param.EVALUATE;
  Cluster[0].radius:= 0;
  Cluster[0].NumSubClusters:= 0;
end;
//==============================================================================
//==============================================================================

procedure TCurve.GENERATE_BASE_CLUSTERS(const safeStep: double);
var

  k, kk: integer;
  NumBaseClusters: integer;
  value, ratio, rat: double;
  max, min: array[0..2] of double;
  Net: TNet1D;
  arcLengthCorrection: double;
begin

  ratio:= 0.5;
  arcLengthCorrection:= 2*arcsin(0.5*ratio)/ratio;
  if length(self.Inward) > 0 then begin

    if self.Inward[0] = true then begin

      self.startPoint.t:= self.ParaPoint[0];
      self.endPoint.t:= self.ParaPoint[1];
    end else begin

      self.startPoint.t:= self.ParaPoint[1];
      self.endPoint.t:= self.ParaPoint[0];
    end;
  end else begin

    self.startPoint.t:= 0;
    self.endPoint.t:= 2*pi;
  end;
  self.startPoint:= self.Param.EVALUATE_FUNCTION_AND_DERIVATIVES(self.startPoint.t);
  self.endPoint:= self.Param.EVALUATE_FUNCTION_AND_DERIVATIVES(self.endPoint.t);

  Net:= MESH_CURVE(self.Param.EVALUATE_FUNCTION_AND_DERIVATIVES, self.startPoint, self.endPoint,
                   self.Param.ParamType = ptStraightLine,
                   safeStep, ratio, ctParametric, true);
  isLoop:= (NORM(DIF(Net[0].point.r, Net[length(Net)-1].point.r)) < 1e-10) and
           (NORM(DIF(Net[0].point.rt, Net[length(Net)-1].point.rt)) < 1e-10);
  if isLoop then begin

    setlength(Net, length(Net)-1);
  end;

  NumBaseClusters:= length(Net);
  setlength(Cluster, NumBaseClusters+1);
  for k:= 0 to NumBaseClusters-1 do begin

    Cluster[k].clusterType:= ctTube;
    Cluster[k].center:= Net[k].r;
    value:= 1/NORM(Net[k].rt);
    Cluster[k].n[0]:= value*Net[k].rt[0];
    Cluster[k].n[1]:= value*Net[k].rt[1];
    Cluster[k].n[2]:= value*Net[k].rt[2];
    // highest distance between adjacent nodes (value)
    if k > 0 then begin

      value:= NORM(DIF(Net[k].r, Net[k-1].r));
      if k < NumBaseClusters-1 then begin

        if NORM(DIF(Net[k+1].r, Net[k].r)) > value then begin

          value:= NORM(DIF(Net[k+1].r, Net[k].r));
        end;
      end;
    end else begin

      value:= NORM(DIF(Net[k+1].r, Net[k].r));
    end;
    rat:= arcLengthCorrection*value*NORM(Net[k].point.rtt)/DOT(Net[k].point.rt, Net[k].point.rt);
    Cluster[k].toln:= 0.5*rat;
    Cluster[k].tolr:= 0.125*value*rat;
    Cluster[k].radius:= 0.55*value;
    Cluster[k].NumSubClusters:= 0;
    if k = 0 then begin

      for kk:= 0 to 2 do begin

        max[kk]:= Net[k].r[kk]+Cluster[k].tolr;
        min[kk]:= Net[k].r[kk]-Cluster[k].tolr;
      end;
    end;
    if k < NumBaseClusters-1 then begin

      for kk:=0 to 2 do begin

        value:= Net[k+1].r[kk];
        if value+Cluster[k+1].tolr > max[kk] then max[kk]:= value+Cluster[k+1].tolr;
        if value-Cluster[k+1].tolr < min[kk] then min[kk]:= value-Cluster[k+1].tolr;
      end;
    end;
  end;
  Cluster[NumBaseClusters].clusterType:= ctSphere;
  Cluster[NumBaseClusters].radius:= 0.5*sqrt(intpower(max[0]-min[0],2)+
                                             intpower(max[1]-min[1],2)+
                                             intpower(max[2]-min[2],2));
  for kk:= 0 to 2 do begin

    Cluster[NumBaseClusters].center[kk]:= 0.5*(min[kk]+max[kk]);
  end;
  Cluster[NumBaseClusters].NumSubClusters:= NumBaseClusters;
  Cluster[NumBaseClusters].FirstSubClusterIndex:= 0;
  // convert TNet1D to TNet
  setlength(self.net, length(Net));
  for k:= 0 to length(Net)-1 do begin

    setlength(self.net[k].p, 2);
    self.net[k].p[0]:= Net[k].point.t;
    self.net[k].p[1]:= sqrt(DOT(Net[k].point.rt, Net[k].point.rt));
    self.net[k].cluster:= Cluster[k];

    setlength(self.net[k].neighbor, 2);
    self.net[k].neighbor[0]:= k-1;
    self.net[k].neighbor[1]:= k+1;
    if isLoop and (k = 0) then begin

      self.net[0].neighbor[0]:= length(Net)-1;
    end;
    if k = length(Net)-1 then begin

      self.net[k].neighbor[1]:= -1;
      if isLoop then begin

        self.net[k].neighbor[1]:= 0;
      end;
    end;
  end;
end;
//==============================================================================

function TCurve.IS_INSIDE(const t: Double; const tol: double): boolean;
var

  ti, tf: double;
  Param: TParam1D;
  p: pointer;
begin

  if length(self.ParaPoint) = 0 then begin

    result:= true;
  end else begin

    p:= self.Param;
    Param:= p;
    if self.Inward[0] = True then begin

      ti:= self.ParaPoint[0];
      tf:= self.ParaPoint[1];
      if ti < tf then begin

        // no parametric jump
        if ((ti-t)*NORM(Param.rt(ti)) < tol) and
           ((t-tf)*NORM(Param.rt(tf)) < tol) then begin

          result:= true;
        end else begin

          result:= false;
        end;
      end else begin

        // parametric jump (e.g., arc)

      end;
    end else begin


    end;
  end;
end;

//==============================================================================

end.
