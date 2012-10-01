

\section{Covariance Model Comparison}

\begin{code}
{-# LANGUAGE NamedFieldPuns, RecordWildCards, PatternGuards,DeriveDataTypeable #-}
module Main where

import Control.Monad
import Data.Array ((!),bounds)
import Data.Function
import Data.List
import Data.Maybe
import Data.Ord
import Debug.Trace
import System.Directory
import System.Environment (getArgs)
import Text.Printf
import Text.Regex.Posix
import System.Console.CmdArgs
import Control.Parallel.Strategies

import qualified Data.MemoCombinators as Memo

import Biobase.Infernal.CM
import Biobase.Infernal.CM.Import
\end{code}



Used for CmdArgs interpretation.

\begin{code}
data Options = Options {global :: Bool, files :: [String]}
  deriving (Show,Data,Typeable)

options = mode $ Options
  { global = False &= text "use global CM comparison"
  , files = def &= args & typ "2 FILES"
  }
\end{code}



Distance measures between two Covariance Models.

\begin{code}
type Distance a =
  ( a
  , a -> a -> a
  , Maybe Double -> a -> a
  , Maybe Double -> a -> a
  , ((String,Double),(String,Double)) -> a -> a
  , ((String,Double),(String,Double)) -> a -> a
  , ((String,Double),(String,Double)) -> a -> a
  , [a] -> [a]
  )

mindist :: Distance Double
mindist = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = 0
  branch = (+)
  transl = trans
  transr = trans
  trans Nothing c = c
  trans (Just t) c = c+t
  mp (s1,s2) c = c + snd s1 + snd s2
  xl (s1,s2) c = c + snd s1 + snd s2
  xr (s1,s2) c = c + snd s1 + snd s2
  optimum [] = []
  optimum xs = [maximum xs]

avgdist :: Distance Double
avgdist = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = 0
  branch = (+)
  transl = trans
  transr = trans
  trans Nothing c = c
  trans (Just t) c = c+t
  mp (s1,s2) c = c + snd s1 + snd s2
  xl (s1,s2) c = c + snd s1 + snd s2
  xr (s1,s2) c = c + snd s1 + snd s2
  optimum []  = []
  optimum xs' =
    let
      --xs = filter (> (-1000000)) xs'
      xs = filter (>= 0) xs'
    in case xs of
      [] -> []
      _  -> [sum xs / genericLength xs]

\end{code}



Parameterized Distance measures.

\begin{code}

--should still use only one trans function ?!
--make sure xs contains no -infty on call
--TODO better yet transform -infty in something like -100000
--TODO EVEN BETTER disallow going into -infty states!!!
--TODO put this idea into the CM import library (with a switch!)
skewL :: Distance (Double,Double)
skewL = pairWith f where
  f xs = [maximumBy (\(a,b) (c,d) -> a `compare` c) xs]

smalldif :: Distance (Double,Double)
smalldif = pairWith f where
  f xs = [maximumBy (\(a,b) (c,d) -> (a+b-(abs $ a-b)) `compare` (c+d-(abs $ c-d))) xs]

maximin :: Distance (Double,Double)
maximin = pairWith f where
  f xs = [maximumBy (\(a,b) (c,d) -> (min a b) `compare` (min c d)) xs]

-- pairs with different optimization function
pairWith f = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = (0,0)
  branch (l1,l2) (r1,r2) = (l1+r1,l2+r2)
  transl Nothing c = c
  transl (Just t) (c1,c2) = (c1+t,c2)
  transr Nothing c = c
  transr (Just t) (c1,c2) = (c1,c2+t)
  mp (s1,s2) (c1,c2) = (c1 + snd s1,c2 + snd s2)
  xl (s1,s2) (c1,c2) = (c1 + snd s1,c2 + snd s2)
  xr (s1,s2) (c1,c2) = (c1 + snd s1,c2 + snd s2)
  optimum [] = []
  optimum xs = f xs

\end{code}



String distances, used in combination with the algebra product operation to
output traces.

\begin{code}

rnaseq :: Distance String
rnaseq = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = ""
  branch = (++)
  transl = trans
  transr = trans
  trans _ c = c
  mp (([b1,b2],_),_) c = [b1] ++ c ++ [b2]
  xl ((b,_),_) c = b ++ c
  xr ((b,_),_) c = c ++ b
  optimum = take 1

bracket :: Distance String
bracket = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = "-"
  branch l r = l ++ "+" ++ r
  transl = trans
  transr = trans
  trans _ c = c
  mp (([b1,b2],_),_) c = "(" ++ c ++ ")"
  xl ((b,_),_) c = "p" ++ c
  xr ((b,_),_) c = c ++ "q"
  optimum = take 1

rnadebug :: Distance String
rnadebug = (init,branch,transl,transr,mp,xl,xr,optimum) where
  init = "-"
  branch l r = l ++ "+" ++ r
  transl = trans
  transr = trans
  trans _ c = c
  mp (([b1,b2],_),_) c = [b1] ++ c ++ [b2]
  xl ((b,_),_) c = b ++ c
  xr ((b,_),_) c = c ++ b
  optimum = take 1

dist1 *** dist2 = (init,branch,transl,transr,mp,xl,xr,optimum) where
  (init1,branch1,transl1,transr1,mp1,xl1,xr1,optimum1) = dist1
  (init2,branch2,transl2,transr2,mp2,xl2,xr2,optimum2) = dist2
  init = (init1,init2)
  branch (a1,a2) (b1,b2) = (branch1 a1 b1, branch2 a2 b2)
  transl t (c1,c2) = (transl1 t c1, transl2 t c2)
  transr t (c1,c2) = (transr1 t c1, transr2 t c2)
  mp (s1,s2) (c1,c2) = (mp1 (s1,s2) c1, mp2 (s1,s2) c2)
  xl (s1,s2) (c1,c2) = (xl1 (s1,s2) c1, xl2 (s1,s2) c2)
  xr (s1,s2) (c1,c2) = (xr1 (s1,s2) c1, xr2 (s1,s2) c2)
  optimum xs = [(x1,x2) | x1 <- nub $ optimum1 [y1 | (y1,y2) <- xs]
                        , x2 <-       optimum2 [y2 | (y1,y2) <- xs, y1 == x1]
               ]
\end{code}



Assymetric state assignment is not easy to score, consider:

S     S
MATP
MATP
MATP
~~~~  ~~~~

any inserted pairs on the left can be "eaten" on the right by IL/IR but we can
not match one state on the left to two states on the right. This would require
going back to nodes...  This means that it is much easier and better to not
consider those states at all.

\begin{code}
recurse (init,branch,transl,transr,mp,xl,xr,optimum) local rfl rfr =
  optimum . concat $ (rec 0 0) : [localbegin i j | local, i <- [0..lastl], j <- [0..lastr], not $ i==0 && j==0] where
  localbegin i j = do
    res <- rec i j
    return . transr localbeginTl $ transl localbeginTr res
  rec = Memo.memo2 (Memo.arrayRange (0,lastl)) (Memo.arrayRange (0,lastr)) rec'
  --the recursive calculation, disallow loops into selfstates
  --[init | local] allows local start anywhere (current (l,r) is empy init then!)
  --rec' :: Int -> Int -> Score
  rec' l r

    | tl=="E" && tr=="E" = [init]

    | tl=="S" && tr=="S" || tl=="D" && tr=="D" = opt $ do
      (l',lT) <- csl ++ localendl
      (r',rT) <- csr ++ localendr
      guard $ l'/=l || r'/=r
      res <- rec l' r'
      return . transr rT $ transl lT res

    | tl=="B" && tr=="B" = opt $
      [ branch rl rr | rl <- rec (left csl)  (left csr),rr <- rec (right csl) (right csr) ] ++
      [ ll `branch` mm `branch` rr | ll <- rec lastl (left csr), mm <- rec (left csl)  (right csr), rr <- rec (right csl) lastr ] ++
      [ ll `branch` mm `branch` rr | mm <- rec (right csl) (left csr), ll <- rec (left csl)  lastr, rr <- rec lastl (right csr) ]

    | tl=="B" = opt $
      [ ll `branch` rr | ll <- rec (left csl) r    , rr <- rec (right csl) lastr ] ++
      [ ll `branch` rr | ll <- rec (left csl) lastr, rr <- rec (right csl) r     ]

    | tr=="B" = opt $
      [ ll `branch` rr | ll <- rec l (left csr)    , rr <- rec lastl (right csr) ] ++
      [ ll `branch` rr | ll <- rec lastl (left csr), rr <- rec l (right csr)     ]

    | tl=="S" = opt $ do
      (l',lT) <- csl
      res <- rec l' r
      return $ transl lT res

    | tr=="S" = opt $ do
      (r',rT) <- csr
      res <- rec l r'
      return $ transr rT res

    | tl `elem` lstates && tr `elem` lstates = opt $ symmetricEmission xl

    | tl `elem` rstates && tr `elem` rstates = opt $ symmetricEmission xr

    | tl=="MP" && tr=="MP" = opt $ symmetricEmission mp

    | tl=="D" && tr `elem` ["E","S"] = do
      (l',lT) <- csl ++ localendl
      res <- rec l' r
      return $ transl lT res

    | tl `elem` ["E","S"] && tr=="D" = do
      (r',rT) <- csr ++ localendr
      res <- rec l r'
      return $ transr rT res

    | otherwise = [] -- trace (show (tl,tr)) []
  {-
    | tl == "B" = opt $
      [ ll `branch` rr | ll <- rec (left csl)  r, rr <- rec (right csl) lastr ] ++
      [ ll `branch` rr | rr <- rec (right csl) r, ll <- rec (left csl)  lastr ]
    | tr == "B" = opt $
      [ ll `branch` rr | ll <- rec l (left csr), rr <- rec lastl (right csr) ] ++
      [ ll `branch` rr | rr <- rec l (right csr), ll <- rec lastl (left csr) ]
    -}
    where
      tl = stype $ ssl ! l
      tr = stype $ ssr ! r
      csl = schildren $ ssl ! l
      csr = schildren $ ssr ! r
      esl = semissions $ ssl ! l
      esr = semissions $ ssr ! r
      sidl = sid $ ssl ! l
      sidr = sid $ ssr ! r
      lstates = ["ML","IL"]
      rstates = ["MR","IR"]
      left xs = fst $ xs !! 0
      right xs = fst $ xs !! 1

      -- some functions
      opt xs = optimum $ parMap rdeepseq id $ [init | local] ++ xs
      symmetricEmission f = do
        (l',lT) <- children rfl l ++ localendl
        (r',rT) <- children rfr r ++ localendr
        guard $ l' /= l || r' /= r
        res <- rec l' r'
        e <- zip esl esr
        return . f e . transr rT $ transl lT res
      -- either l or r does an emission
      asymmetricEmission f = do
        (l',lT) <- children rfl l ++ localendl
        (r',rT) <- children rfr r ++ localendr
        guard $ l' /= l || r' /= r
        res <- rec l' r'
        let l = max (length esl) (length esr)
        let pbs = map fst $ if length esl == l then esl else esr
        e <- zip (zip pbs (map snd esl ++ repeat penalty)) (zip pbs (map snd esr ++ repeat penalty))
        return . f e . transr rT $ transl lT res

  -- important stuff
  ssl = states rfl
  ssr = states rfr
  lastl = snd $ bounds ssl
  lastr = snd $ bounds ssr
  nodesL = fromIntegral . snd . bounds $ nodes rfl
  nodesR = fromIntegral . snd . bounds $ nodes rfr
  -- how much to penalize asymmetric insertion of nucleotides
  penalty = -1
  -- allow cheap transition to "E"nd state
  -- we cheat a bit by making the transition slightly cheaper (dividing by nodes, not states)
  -- this makes sense as transitions into anything but MP,ML,MR should give bad scores anyway
  localendl = [(lastl,Just $ log (default_pend / nodesL) / log 2) | local]
  localendr = [(lastr,Just $ log (default_pend / nodesR) / log 2) | local]
  -- transition score for local beginnings
  localbeginTl = Just $ log (default_pbegin / nodesL) / log 2
  localbeginTr = Just $ log (default_pbegin / nodesR) / log 2
  -- default values to use for local begin/end states
  default_pbegin = 0.05
  default_pend   = 0.05

childrenOnly rf k = map fst . schildren $ states rf !  k
\end{code}



This one is tricky: always return the children, except if we are an end state,
then we return ourselves to allow global alignment.

integrating "local" as a direct jump to E with no transition cost should be
better (no trans cost into the local start state)

\begin{code}
children rf k =
  let
    state = states rf ! k
  in
    if stype state == "E" && False
    then [(k,Nothing)]
    else schildren state
\end{code}



\begin{code}
mean s1 s2 sl sr = (s1+s2) / (sl + sr)
tmax s1 s2 sl sr = max (s1/sl) (s2/sr)

main = do
  -- (a:b:[]) <- getArgs
  ca <- cmdArgs "CMCompare (c) Christian Hoener zu Siederdissen, 2009, choener@tbi.univie.ac.at" [options]
  isloud <- isLoud
  isquiet <- isQuiet
  isnormal <- isNormal
  if (length (files ca) /= 2)
    then do
      print "give exactly TWO CMs, ... please!"
    else do
      let a = files ca !! 0
      let b = files ca !! 1
      ae <- doesFileExist a
      be <- doesFileExist b
      if (ae && be)
        then do
        Right rf1' <- fromFile a
        let rf1 = canonize rf1'
        Right rf2' <- fromFile b
        let rf2 = canonize rf2'

        let both = recurse (maximin *** (rnaseq *** (rnadebug *** bracket))) (not $ global ca) rf1 rf2
        let jrf1 = map (/2) $ recurse mindist False rf1 rf1
        let jrf2 = map (/2) $ recurse mindist False rf2 rf2

        if (null both || null jrf1 || null jrf2)
          then do
          putStrLn . concat $ intersperse " " ["error with: ", a, b, show both, show jrf1, show jrf2]
          else do
          let ((sl,sr),(rnaseq,(rnadebug,bracket))) = head both
          let jl = head jrf1
          let jr = head jrf2
          let res = tmax sl sr jl jr
          -- printf "Identity: %7.2f %7.2f %7.2f %7.2f %7.2f %s %s   %s %s\n" res sl sr jl jr rna bracket a b
          when (isquiet && not isnormal && not isloud) (printf "Comparison: %7.2f %s %s\n" res a b)
          when (isnormal) (printf "Comparison: %7.2f %7.2f %7.2f %7.2f %7.2f %s   %s %s\n" res sl sr jl jr rnaseq a b)
          when (isloud) $ do
            printf "Verbose: %s\n" rnadebug
            printf "Verbose: %s\n" bracket
        else do
          when (not ae) (putStrLn $ "missing file: " ++ a)
          when (not be) (putStrLn $ "missing file: " ++ b)
\end{code}
