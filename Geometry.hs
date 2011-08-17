{-# LANGUAGE TypeSynonymInstances #-}

module Geometry where

import Data.List
import Data.Maybe

import Control.Monad

data Vector = Vector !Double !Double !Double
  deriving (Show, Read)

instance Eq Vector where
  Vector x1 y1 z1 == Vector x2 y2 z2 = x1 ~= x2 && y1 ~= y2 && z1 ~= z2

(~=) :: (Floating a, Ord a) => a -> a -> Bool
a ~= b = let c = abs (a - b) 
         in -0.00001 < c && c < 0.00001

(~~) :: (Floating a, Ord a) => a -> a -> Bool
a ~~ b | a == 0 = -0.001 < b && b < 0.001
       | b == 0 = -0.001 < a && a < 0.001
       | otherwise = let c = a / b 
                     in 0.999 < c && c < 1.001

instance Num Vector where
  Vector x1 y1 z1 + Vector x2 y2 z2 = Vector (x1 + x2) (y1 + y2) (z1 + z2)
  Vector x1 y1 z1 - Vector x2 y2 z2 = Vector (x1 - x2) (y1 - y2) (z1 - z2)
  Vector x1 y1 z1 * Vector x2 y2 z2 = Vector (y1 * z2 - z1 * y2) (z1 * x2 - x1 * z2) (x1 * y2 - y1 * x2)
  negate (Vector x y z) = Vector (-x) (-y) (-z)
  abs    (Vector x y z) = Vector (abs x) (abs y) (abs z)
  signum v = (1 / norm v) |*| v
  fromInteger n = Vector (fromInteger n) (fromInteger n) (fromInteger n)

unit :: Vector -> Vector
unit = signum

(|*|) :: Double -> Vector -> Vector
k |*| Vector x y z = Vector (k * x) (k * y) (k * z)

norm :: Vector -> Double
norm (Vector x y z) = sqrt (x^2 + y^2 + z^2)

scaleTo :: Vector -> Double -> Vector
v `scaleTo` l = l |*| unit v

distance :: Point -> Point -> Double
distance p1 p2 = norm (p1 - p2)

dot :: Vector -> Vector -> Double
Vector x1 y1 z1 `dot` Vector x2 y2 z2 = x1 * x2 + y1 * y2 + z1 * z2

type Angle = Double

angleBetween :: Vector -> Vector -> Angle
angleBetween u v = acos (clamp (-1,1) (unit u `dot` unit v))

basicAngle :: Angle -> Angle
basicAngle ang | ang < 0   = basicAngle (ang + tpi)
               | ang > tpi = basicAngle (ang - tpi)
               | otherwise = ang
  where tpi = 2 * pi

type Point = Vector

origin :: Point
origin = Vector 0 0 0

type Polygon = [Point]

midpoint :: Point -> Point -> Point
midpoint a b = 0.5 |*| (a + b)

centroid :: Polygon -> Point
centroid [] = origin
centroid ps = (1 / fromIntegral (length ps)) |*| sum ps

polygonNormal :: Polygon -> Vector
polygonNormal (a:b:c:_) = unit ((b - a) * (c - a))


class Transform3D a where
  translate3D :: Vector          -> a -> a
  rotate3D    :: Vector -> Angle -> a -> a
  scale3D     :: Double          -> a -> a

instance Transform3D Vector where
  translate3D trans = (+ trans)
  rotate3D z t v = (cos t |*| v) + (sin t |*| (z * v)) + (((z `dot` v) * (1 - cos t)) |*| z) -- Rodrigue's Formula
  scale3D = (|*|)

instance (Transform3D a) => Transform3D [a] where
  translate3D v     = map (translate3D v)
  rotate3D axis ang = map (rotate3D axis ang)
  scale3D k         = map (scale3D k)

instance (Transform3D a, Transform3D b) => Transform3D (a,b) where
  translate3D v     (a,b) = (translate3D v     a, translate3D v     b)
  rotate3D axis ang (a,b) = (rotate3D axis ang a, rotate3D axis ang b)
  scale3D k         (a,b) = (scale3D k         a, scale3D k         b)


alignToTwist :: (Transform3D a) => Vector -> Vector -> Angle -> a -> a
alignToTwist v' v1 twist a = rotate3D (rotAxis v1 v') (angleBetween v1 v' + twist) a

alignTo :: (Transform3D a) => Vector -> Vector -> a -> a
alignTo v' v1 a = alignToTwist v' v1 0 a

rotAxis v1 v' = let axis' = v1 * v'
                in if axis' == Vector 0 0 0
                    then let Vector vx vy vz = v1
                             vec | vx == 0   = Vector 0 vz (-vy)
                                 | vy == 0   = Vector vz 0 (-vx)
                                 | otherwise = Vector vy (-vx) 0
                         in unit vec
                    else unit axis'

type Rotation = (Vector, Angle)

identity :: Rotation
identity = (Vector 0 0 1, 0)


class (Ord a) => Interval a where
  clamp :: (a,a) -> a -> a
  transformToInterval :: (a,a) -> (a,a) -> a -> a  

instance Interval Double where
  clamp (a,b) c | c < a = a
                | c > b = b
                | otherwise = c
  transformToInterval (a,b) (c,d) e = let dInt = d - c
                                      in if dInt == 0
                                          then b
                                          else a + (e - c) * (b - a) / dInt
instance Interval Float where
  clamp (a,b) c | c < a = a
                | c > b = b
                | otherwise = c
  transformToInterval (a,b) (c,d) e = let dInt = d - c
                                      in if dInt == 0
                                          then b
                                          else a + (e - c) * (b - a) / dInt

data HullPoint = HullPoint {hullPosition :: Point, hullRadius :: Double}
  deriving (Show, Eq, Read)

instance Transform3D HullPoint where
  translate3D v     (HullPoint p r) = HullPoint (translate3D v     p) r
  rotate3D axis ang (HullPoint p r) = HullPoint (rotate3D axis ang p) r
  scale3D k         (HullPoint p r) = HullPoint (scale3D k         p) r


type HullSegment = (HullPoint, HullPoint)
type HullPolygon = [HullPoint]


data Plane = Plane Vector Double
  deriving (Show, Eq)

instance Transform3D Plane where
  translate3D v     (Plane n c) = Plane n (c + v `dot` n)
  rotate3D axis ang (Plane n c) = Plane (rotate3D axis ang n) c
  scale3D k         (Plane n c) = if k < 0
                                   then Plane (-n) (-c * k)
                                   else Plane n (c * k)


data HullFace = HullFace {hullFacePolygon :: HullPolygon, hullFacePlane :: Plane}
  deriving (Show, Eq) 

instance Transform3D HullFace where
  translate3D v     (HullFace hps pl) = HullFace (translate3D v     hps) (translate3D v     pl)
  rotate3D axis ang (HullFace hps pl) = HullFace (rotate3D axis ang hps) (rotate3D axis ang pl)
  scale3D k         (HullFace hps pl) = HullFace (scale3D k         hps) (scale3D k         pl)


data ConvexHull = ConvexHull [HullPoint] [HullSegment] [HullFace]
  deriving (Show, Eq)

instance Transform3D ConvexHull where
  translate3D v     (ConvexHull hps hss hfs) = ConvexHull (translate3D v     hps) (translate3D v     hss) (translate3D v     hfs)
  rotate3D axis ang (ConvexHull hps hss hfs) = ConvexHull (rotate3D axis ang hps) (rotate3D axis ang hss) (rotate3D axis ang hfs)
  scale3D k         (ConvexHull hps hss hfs) = ConvexHull (scale3D k         hps) (scale3D k         hss) (scale3D k         hfs)


atan3 y x = if x == 0
             then if y > 0
                   then pi/2
                   else 3*pi/2
             else if x > 0
                   then if y > 0
                         then atan (y/x)
                         else atan (y/x) + 2*pi
                   else atan (y/x) + pi

threeSpheresTangentPlane :: HullPoint -> HullPoint -> HullPoint -> Vector -> Plane
threeSpheresTangentPlane (HullPoint p1' r1') (HullPoint p2' r2') (HullPoint p3' r3') v =
  if r1' ~= r2' && r2' ~= r3'
   then let zv = unit ((p2' - p1') * (p3' - p1'))
            n  = if zv `dot` v > 0
                  then  zv
                  else -zv
        in Plane n ((n `dot` p1') + r1')
   else let [(p1,r1),(p2,r2),(p3,r3)] | r1' ~= r2' = [(p3',r3'),(p1',r1'),(p2',r2')]
                                      | r2' ~= r3' = [(p1',r1'),(p2',r2'),(p3',r3')]
                                      | r3' ~= r1' = [(p2',r2'),(p3',r3'),(p1',r1')]
                                      | otherwise  = sortBy (\(_,a) (_,b) -> compare (-a) (-b)) [(p1',r1'),(p2',r2'),(p3',r3')]
            sc2 = similarityCenter p1 r1 p2 r2
            sc3 = similarityCenter p1 r1 p3 r3
            d1  = norm ((sc2 - p1) * (sc3 - p1)) / distance sc2 sc3
            ang = asin (r1 / d1)
            zv  = unit ((p2 - p1) * (p3 - p1))
            n   = if zv `dot` v > 0
                   then rotate3D (unit (sc3 - sc2)) ang   zv
                   else rotate3D (unit (sc2 - sc3)) ang (-zv)
        in Plane n (n `dot` sc2)

similarityCenter p1 r1 p2 r2 =
  let k = distance p1 p2
      l = r1 * k / (r2 - r1)
  in p2 + ((p1 - p2) `scaleTo` (l + k))

insideHalfSpace :: HullPoint -> HullFace -> Bool
HullPoint p r `insideHalfSpace` HullFace _ (Plane v k) =
  let a = p `dot` v
      b = k - r
  in a ~= b || a < b

inConeOf :: HullPoint -> HullPoint -> HullPoint -> Bool
inConeOf (HullPoint p2' r2') (HullPoint p1' r1') (HullPoint p r) =
  if r1' ~~ r2'
   then d + r <= r1'
   else let (p2,r2,p1,r1) = if r1' < r2'
                             then (p2',r2',p1',r1')
                             else (p1',r1',p2',r2')
            sc = similarityCenter p2 r2 p1 r1
            coneSine = r1 / distance p1 sc
            sc' = sc + ((p1 - sc) `scaleTo` (r / coneSine))
            testSine = (norm ((p2 - p) * (p1 - p)) / (distance p2 p1 * distance p sc')) 
        in testSine <= coneSine
  where d = norm ((p2' - p) * (p1' - p)) / (distance p2' p1')

coplanarHullPoints :: HullPoint -> HullPoint -> Plane -> [HullPoint] -> Bool -> (HullFace, [HullFace])
coplanarHullPoints hp1@(HullPoint p1 r1) hp2 pl@(Plane n _) hps starter =
  let hpAll = hp1:hp2:hps
      cen = centroid (map hullPosition hpAll)
      v = hullPosition hp2 - p1
      u = v * (v - n)
      (an, Plane sn _) = if u `dot` (p1 - cen) > 0
                          then let an' = n
                               in (an', threeSpheresTangentPlane hp1 (HullPoint (p1 + an') r1) hp2   u ) 
                          else let an' = -n
                               in (an', threeSpheresTangentPlane hp1 (HullPoint (p1 + an') r1) hp2 (-u))
      sorted@(x:xs) = wrap2D hp1 cen an hp1 hp2 sn hpAll
      hfs = if starter
             then getHullFaces (sorted ++ [x, head xs])
             else removeEdge (hp1,hp2) (getHullFaces (sorted ++ [x, head xs]))
  in (HullFace sorted pl, hfs)
  where sortRing v@(Vector x y z) hps =
          let i = unit (v * Vector x y (2 * maximum [x,y,z]))
              j = v * i
              pairs = map (projectPairup i j) hps
          in map snd (sortBy (\(a,_) (b,_) -> compare a b) pairs)
          where projectPairup i j hp@(HullPoint p _) = (atan3 (p `dot` j) (p `dot` i), hp)

        getHullFaces (a:xs@(b:c:_)) = HullFace [a,b,c] pl : getHullFaces xs
        getHullFaces _ = []

        removeEdge (a,b) = removeBy ((\x -> x == [a,b] || x == [b,a]) . (take 2) . hullFacePolygon)

        wrap2D hp0 cen an hp1 hp2@(HullPoint p2 r2) sn hps =
          let p2' = p2 + an
              pairs = mapMaybe (pairWithAng p2 p2' sn hp2 (HullPoint p2' r2)) (hps \\ [hp1,hp2])
              (_,(n,hp)) = minimumBy (\(a,_) (b,_) -> compare a b) pairs
          in if hp == hp0
              then [hp1,hp2]
              else  hp1 : wrap2D hp0 cen an hp2 hp n hps
          where pairWithAng p1 p2 v hp1 hp2 hp@(HullPoint p _) =
                  if inConeOf hp1 hp2 hp
                   then Nothing
                   else let Plane nv _ = threeSpheresTangentPlane hp1 hp2 hp ((p1 - p) * (p2 - p))
                        in Just (angleBetween v nv, (nv,hp))

removeBy f (x:xs) = if f x
                     then xs
                     else x : removeBy f xs
removeBy _ _ = []

convexHull :: [HullPoint] -> ConvexHull
convexHull [hp] = ConvexHull [hp] [] []
convexHull [hp1,hp2] = ConvexHull [hp1,hp2] [(hp1,hp2)] []
convexHull hps' =
  let hps = nub hps'
      xPairs = map pairWithXR hps
      (_,hp1) = minimumBy (\(a,_) (b,_) -> compare a b) xPairs

      HullPoint p1@(Vector x y z) r1 = hp1
      tp2 = Vector x y (z + 2*r1)
      hps1 = delete hp1 hps
      angPairs = mapMaybe (pairWithAng p1 tp2 (Vector (-1) 0 0) hp1 (HullPoint tp2 r1)) hps1
      (_,(nv,hp2)) = minimumBy (\(a,_) (b,_) -> compare a b) angPairs

      p2 = hullPosition hp2
      hps2 = delete hp2 hps1
      plPairs = mapMaybe (pairWithPlanes p2 p1 nv hp2 hp1) hps2
      sorted = sortBy (\(a,_) (b,_) -> compare a b) plPairs
      (ang,(pl,_)) = head sorted
      cohps = map (snd . snd) (takeWhile ((~= ang) . fst) sorted)
      (hf0,hfs) = coplanarHullPoints hp1 hp2 pl cohps True
  in giftWrapping (hullFacePolygon hf0) [] [hf0] hfs hps
  where pairWithXR hp@(HullPoint (Vector x _ _) r) = (x - r, hp)

        pairWithAng p1 p2 xv hp1 hp2 hp@(HullPoint p _) =
          if inConeOf hp1 hp2 hp
           then Nothing
           else let pl@(Plane nv _) = threeSpheresTangentPlane hp1 hp2 hp ((p - p2) * (p - p1))
                in Just (angleBetween xv nv, (nv,hp))

        pairWithPlanes p1 p2 v hp1 hp2 hp@(HullPoint p _) =
          if inConeOf hp1 hp2 hp
           then Nothing
           else let pl@(Plane nv _) = threeSpheresTangentPlane hp1 hp2 hp ((p - p2) * (p - p1))
                in Just (angleBetween v nv, (pl,hp))

giftWrapping :: [HullPoint] -> [HullSegment] -> [HullFace] -> [HullFace] -> [HullPoint] -> ConvexHull
giftWrapping hps hss hfs (HullFace [hp1,hp2,hp3] (Plane n _) : hsfs) points =
  let HullPoint p1 r1 = hp1
      HullPoint p2 r2 = hp2
      p3 = hullPosition hp3
      tn = (p1 - p3) * (p2 - p3)
      k = unit (p2 + (n `scaleTo` r2) - p1 - (n `scaleTo` r1))
      rest = points \\ [hp1,hp2]
      pairs = if tn `dot` n > 0
               then mapMaybe (pairWithPlanes p2 p1 n ( k * n)) rest
               else mapMaybe (pairWithPlanes p1 p2 n (-k * n)) rest
      sorted = sortBy (\(a,_) (b,_) -> compare a b) pairs
      (ang,(pl@(Plane nn _),ppp)) = head sorted
      cohps = map (snd . snd) (takeWhile ((~= ang) . fst) sorted)
      (hf1,hfs1) = coplanarHullPoints hp1 hp2 pl cohps False
  in if any (doneBefore (hp1,hp2) nn) hsfs
      then giftWrapping hps ((hp1,hp2):hss) hfs (removeBy (doneBefore (hp1,hp2) nn) hsfs) points
      else giftWrapping (hullFacePolygon hf1 ++ hps) ((hp1,hp2):hss) (hf1:hfs) (hfs1 ++ hsfs)  points
  where pairWithPlanes p1 p2 i j hp@(HullPoint p _) =
          if inConeOf hp1 hp2 hp
           then Nothing
           else let pl@(Plane nv _) = threeSpheresTangentPlane hp1 hp2 hp ((p1 - p) * (p2 - p))
                    ang' = atan3 (nv `dot` j) (nv `dot` i)
                    ang = if ang' ~= (2*pi)
                           then 0
                           else ang'
                in Just (ang, (pl,hp))

        getNormal (HullFace _ (Plane n _)) = n

        doneBefore (a,b) nn (HullFace [c,d,_] (Plane n _)) = (n == nn) && ((a == c && b == d) || (a == d && b == c))

giftWrapping hps hss hfs _ _ = ConvexHull (nub hps) hss hfs