{-# LANGUAGE PatternGuards #-}

module Main where

import Data.List
import Text.Printf
import Data.Tuple.All



-- | Either add a score or write a "<20" to an edge

addScore wp (x,y)
  | Just vs <- (x,y) `lookup` wp = (x,y, (vs !!0) `xmin` (vs!!1))
  | otherwise = (x,y,"<20")

xmin a b = if (read a :: Double) < read b then a else b

-- | ???

outScore wp x =
  let
    hits = map (\((k1,k2),ys) -> (if k1 == x then k2 else k1,(ys !!0) `xmin` (ys!!1))) $ filter (\((k1,k2),_) -> x == k1 || x == k2) wp
    besthit = maximumBy (\(_,x) (_,y) -> (read x :: Double) `compare` read y) hits
  in
    if null hits
    then (x,("all","<20"))
    else (x,besthit)



-- | transform the clan data structure into a dot-readable file

clanGraph (c,ms,outs) = unlines $
  [printf "graph %s {" c] ++
  nodes ++
  outnodes ++
  ["}"]
  where
    f "<20" = "<20"
    f x = printf "%.0f" (read x :: Double)
    nodes = map (\(k1,k2,v) -> printf "%s -- %s [labelfloat=true,label = \"%s\"];" k1 k2 (f v)) ms
    outnodes = map (\(k1,(k2,v)) -> printf "%s -- \"%s\" [label = \"%s\"];" k1 (k1++"\\nvs\\n"++k2) (f v)) outs



-- | prepare a data structure for each clan

clan2graph wp (c,xs) =
  let
    -- add score for each clan edge
    pairs = map (addScore wp) [(x,y) | x <- xs, y <- xs, x > y]
    -- the links to something not in the clan
    outgroup = filter (\((k1,k2),_) -> (k1 `elem` xs && k2 `notElem` xs) || (k1 `notElem` xs && k2 `elem` xs)) wp
    outhits = map (outScore outgroup) xs
  in
    (c,pairs,outhits)



-- |

main = do
  wp <- readFile "rfam-weakpairs"
    >>= return
    . map (\(x:y:xs) -> let x' = take 7 x; y' = take 7 y in if x > y then ((x',y'),xs) else ((y',x'),xs))
    . map words
    . lines
  cnts <- readFile "CLANDESC.all"
    >>= return
    . map (\(x:xs) -> (x,map init xs))
    . groupBy (\x y -> take 2 y /= "CL")
    . map (head . tail)
    . filter (\(x:_) -> x == "MB" || x == "AC")
    . map words
    . lines
  mapM_ (putStrLn . clanGraph . clan2graph wp) cnts
