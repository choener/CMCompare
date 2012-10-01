import Data.List

import System.Environment (getArgs)
import System.Directory
import Text.Regex.Posix
import qualified Data.Set as S

intoParts p [] = []
intoParts p xs = take p xs : intoParts p (drop p xs)

swap (x,y) = if x < y then (x,y) else (y,x)

main = do
  (a:[]) <- getArgs
  dirc <- getDirectoryContents a
  let dc = filter (=~ ".*\\.cm") dirc
  let allall = S.toList . S.fromList $ map swap [ (l,r) | l <- dc, r <- dc ]
  mapM_ (putStrLn . show)  $ intoParts 100 allall
