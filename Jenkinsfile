def artserver = Artifactory.server('store.terradue.com')
def buildInfo = Artifactory.newBuildInfo()
buildInfo.env.capture = true

pipeline {

  options {
    buildDiscarder(logRotator(numToKeepStr: '5'))
  }

  agent { 
    node { 
      label 'ci-community-docker' 
    }
  }

  stages {

    stage('Package & Dockerize') {
      steps {
        
        // See Jenkins's "Global Tool Configuration"
        withMaven( maven: 'apache-maven-3.0.5' ) {
            sh 'mvn -B deploy'
        }

        script {
          artserver.publishBuildInfo buildInfo
        }
      }
    }
  }
}